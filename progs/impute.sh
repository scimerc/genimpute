#!/usr/bin/env bash

set -Eeou pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

declare -r tmpprefix=${opt_outprefix}_tmp

declare -r debuglogfn=${tmpprefix}_debug.log
declare -r sampleprefix=${tmpprefix}_samples/sample
declare -r scriptprefix=${tmpprefix}_scripts/script
declare -r scriptlogprefix=${tmpprefix}_logs/script

mkdir -p $( dirname ${sampleprefix} )
mkdir -p $( dirname ${scriptprefix} )
mkdir -p $( dirname ${scriptlogprefix} )

declare -r cfg_execmode=$( cfgvar_get execmode )
declare -r cfg_igroupsize=$( cfgvar_get igroupsize )
declare -r cfg_genomemap=$( cfgvar_get genomemap )
declare -r cfg_queuecmd=$( cfgvar_get queuecmd )
declare -r cfg_refprefix=$( cfgvar_get refprefix )
declare -r chromosomes="$( seq 1 22 ) 23 25"
declare -r genmap=${BASEDIR}/lib/data/${cfg_genomemap}

bcftoolsexec=${BASEDIR}/lib/3rd/bcftools
eaglexec=${BASEDIR}/lib/3rd/eagle
minimac_version=3
plink_mem=30000

timexec='/usr/bin/time -o /dev/stdout -f "ResStats\tRealSec:\t%e\tMaxMemKB:\t%M\tUserSec:\t%U\tSysSec:\t%S"'

#-------------------------------------------------------------------------------

# purge fam file ids of any unwanted characters

if [ -s ${tmpprefix}_idfix.fam -a -s ${tmpprefix}_idfix.bim -a -s ${tmpprefix}_idfix.bed ] ; then
  echo 'found idfix files. skipping ID purge..'
else
  awk -F $'\t' -f ${BASEDIR}/lib/awk/idclean.awk --source '{
    OFS="\t"
    fid=idclean($1)
    iid=idclean($2)
    gsub( "_+$", "", fid )
    gsub( "^_+", "", iid )
    print( $1, $2, fid, iid )
  }' ${opt_inprefix}.fam > ${tmpprefix}.idmap
  Nold=$( sort -u -k 1,2 ${tmpprefix}.idmap | wc -l )
  Nnew=$( sort -u -k 3,4 ${tmpprefix}.idmap | wc -l )
  if [ $Nold -eq $Nnew ] ; then
    ${plinkexec} --bfile ${opt_inprefix} \
                 --update-ids ${tmpprefix}.idmap \
                 --make-bed --out ${tmpprefix}_idfix
  else
    echo '====> warning: could not purge IDs due to conflicts.'
    cp ${opt_inprefix}.bed ${tmpprefix}_idfix.bed
    cp ${opt_inprefix}.bim ${tmpprefix}_idfix.bim
    cp ${opt_inprefix}.fam ${tmpprefix}_idfix.fam
  fi
fi

#-------------------------------------------------------------------------------

# convert input files to vcf formats

for chr in ${chromosomes} ; do
  if [ -e ${tmpprefix}_chr${chr}.bcf ] ; then
    printf "input file '%s' ready for phasing. skipping conversion..\\n", ${tmpprefix}_chr${chr}.bcf
    continue
  fi
  ${timexec} ${plinkexec} \
    --bfile ${tmpprefix}_idfix \
    --chr ${chr} \
    --recode vcf bgz \
    --out ${tmpprefix}_chr${chr}
  ${bcftoolsexec} convert -Ob --threads 3 \
                  ${tmpprefix}_chr${chr}.vcf.gz \
                  > ${tmpprefix}_chr${chr}_out.bcf
  ${bcftoolsexec} annotate --rename-chrs <( echo -e '23 X\n25 X' ) -Ob --threads 3 \
                  ${tmpprefix}_chr${chr}_out.bcf \
                  > ${tmpprefix}_chr${chr}_outx.bcf
  ${bcftoolsexec} index ${tmpprefix}_chr${chr}_outx.bcf
  mv ${tmpprefix}_chr${chr}_outx.bcf     ${tmpprefix}_chr${chr}.bcf
  mv ${tmpprefix}_chr${chr}_outx.bcf.csi ${tmpprefix}_chr${chr}.bcf.csi
done

#-------------------------------------------------------------------------------

# split individuals in samples

for chr in ${chromosomes} ; do
  ${bcftoolsexec} query -l ${tmpprefix}_chr${chr}.vcf.gz
done | sort -u | split -d -l ${cfg_igroupsize} /dev/stdin ${sampleprefix}
groupsamplesize=${cfg_igroupsize}
tmpsamplelist=$( ls "${sampleprefix}"* )
totalsamplesize=$( cat "${sampleprefix}"* | wc -l )
if [ ${totalsamplesize} -lt ${groupsamplesize} ] ; then
  groupsamplesize=${totalsamplesize}
fi

#-------------------------------------------------------------------------------

# write phase scripts

for chr in ${chromosomes} ; do
chrtag=${chr}
if [ ${chr} -eq 23 -o ${chr} -eq 25 ] ; then
  chrtag=X
fi
#TODO: implement different running time for reference-less phasing (stage ii)
runtimehrs=$(( 6 * ( totalsamplesize / 10000 + 1 ) ))
cat > ${scriptprefix}1_phase_chr${chr}.sh << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G
#SBATCH --time=${runtimehrs}:00:00

# Eagle's parallelization is very efficient. colossus machines seem to have
# 20 cpus but 40 hyperthreading cores. Therefore we use num_cpus*2 threads.
# Running time: 10c-2G-20t From 18min to 3h for batch3 (~9k individuals)

set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_chr${chr}_phased.vcf.gz" ]; then
  printf "phased haplotypes present. nothing to do..\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
${timexec} ${eaglexec} \\
  --chrom ${chrtag} \\
  --geneticMapFile ${genmap} \\
  --vcfRef ${cfg_refprefix}.chr${chr}.haplotypes.bcf \\
  --vcfTarget ${tmpprefix}_chr${chr}.bcf \\
  --outPrefix ${tmpprefix}_chr${chr}_phasing \\
  --numThreads \$(( num_cpus * 2))
mv ${tmpprefix}_chr${chr}_phasing.vcf.gz ${tmpprefix}_chr${chr}_phased.vcf.gz
EOI
done

#-------------------------------------------------------------------------------

# write impute scripts

SBATCH_CONF_MM3="
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=14G
#SBATCH --time=$(( 12 * ( groupsamplesize / 500 + 1 ) )):00:00

# results for batch3; n=9200 individuals total, 14 samples with n=650 each; hrc-32k reference
# Minimac3-omp
#   Number of individuals does not seem to increase mem usage.
#   Mem usage seems to be determined by ref file.
#   3) sbatch: cpus-per-task=4; mem-per-cpu=14 mm3: lowMem; cpus=8  -> 10h; 54jobs; chr2 max.mem=52/56G
#   2) sbatch: cpus-per-task=4; mem-per-cpu=8  mm3: lowMem; cpus=4  -> 14h; 46jobs; chr2 max.mem=31/32G(!)
#   1) sbatch: cpus-per-task=2; mem-per-cpu=16 mm3: lowMem; cpus=2  -> 21h; 46jobs; chr2 max.mem=23/32G
"
SBATCH_CONF_MM4="
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=$(( 4 * ( groupsamplesize / 500 + 1 ) )):00:00

# results for batch3; n=9200 individuals total, 14 samples with n=650 each; hrc-32k reference
# minimac4
#   5) sbatch: cpus-per-task=4 mem-per-cpu=4G  mm4: cpus=8  -> 1.5h; 178jobs; chr6 max.mem=7/16G
#   4) sbatch: cpus-per-task=2 mem-per-cpu=8G  mm4: cpus=4  -> 3.0h; 126jobs; chr6 max.mem=6/16G
"

for chr in ${chromosomes} ; do
  for samplefile in ${tmpsamplelist} ; do
    sampletag=${samplefile##*\/}
    sampletag=${sampletag//*sample/s}
    sampletag=${sampletag//*s/sample}
    imputescriptfn="${scriptprefix}2_impute_chr${chr}_${sampletag}.sh"
    case "${minimac_version}" in
      "3")
        minimacexec="${BASEDIR}/lib/3rd/minimac3_march-sb_omp --lowMemory"
        sbatch_conf=${SBATCH_CONF_MM3}
        ;;
      "4")
        minimacexec=${BASEDIR}/lib/3rd/minimac4
        sbatch_conf=${SBATCH_CONF_MM4}
        ;;
      *)
        printf "error unknown minimac_version '%s'\n" "${minimac_version}" >&2
        exit 1
        ;;
    esac
    cat >> ${imputescriptfn} << EOI
#!/usr/bin/env bash
${sbatch_conf}

set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_chr${chr}_${sampletag}_imputed.dose.vcf.gz" ]; then
  printf "final vcf files present. nothing to do..\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
${bcftoolsexec} view -S ${samplefile} -Oz \\
  --force-samples ${tmpprefix}_chr${chr}_phased.vcf.gz \\
  > ${tmpprefix}_chr${chr}_${sampletag}_phased.vcf.gz
${timexec} ${minimacexec} \\
  --cpus \$(( num_cpus * 2 )) \\
  --haps ${tmpprefix}_chr${chr}_${sampletag}_phased.vcf.gz \\
  --refHaps ${cfg_refprefix}.chr${chr}.m3vcf.gz \\
  --noPhoneHome \\
  --prefix ${tmpprefix}_chr${chr}_${sampletag}_imputing \\
  > ${tmpprefix}_chr${chr}_${sampletag}_imputing.log 2>&1
${bcftoolsexec} index -f ${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.vcf.gz
# convert to bcf
${bcftoolsexec} view ${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.vcf.gz \
    -Ob --threads 3 \
  > ${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.bcf
# fix minimac3 output
${bcftoolsexec} view -h ${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.bcf \\
  | awk -v genofilter='##FILTER=<ID=GENOTYPED,Description="Site was genotyped">' \\
    'BEGIN{ flag = 0 } {
      if ( \$0 ~ "^##FILTER" ) {
        if ( flag==0 ) print( genofilter )
        flag = 1
      }
      print
    }' \\
  > ${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch
${bcftoolsexec} reheader -h ${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch \\
  ${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.bcf \\
  > ${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch.dose.bcf
${bcftoolsexec} index -f ${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch.dose.bcf
mv \\
  ${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch.dose.bcf \\
  ${tmpprefix}_chr${chr}_${sampletag}_imputed.dose.bcf
mv \\
  ${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch.dose.bcf.csi \\
  ${tmpprefix}_chr${chr}_${sampletag}_imputed.dose.bcf.csi
EOI
  done
done

#-------------------------------------------------------------------------------

# write merge scripts

for chr in ${chromosomes} ; do
  cat > ${scriptprefix}3_merge_chr${chr}.sh << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00

set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${opt_outprefix}_chr${chr}_imputed.dose.vcf.gz" ]; then
  printf "merged vcf files present. nothing to do..\\n"
  exit 0
fi
Ns=\$( ls ${tmpprefix}_chr${chr}_*_imputed.dose.bcf | wc -l ) 
if [ \${Ns} -eq 1 ] ; then
  cp ${tmpprefix}_chr${chr}_*_imputed.dose.bcf ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz
else 
  ${bcftoolsexec} merge ${tmpprefix}_chr${chr}_*_imputed.dose.bcf --threads 3 -Oz \\
    > ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz
fi
mv \\
  ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz \\
  ${opt_outprefix}_chr${chr}_imputed.dose.vcf.gz
EOI
done

#-------------------------------------------------------------------------------

# exec

if [ ${opt_dryimpute} -eq 1 ] ; then
  exit 0
fi

echo 'submitting scripts..'

jobscripts=$(ls ${scriptprefix}* | sort -V ) # sort -V places chr1 before chr10

case ${cfg_execmode} in
  "local")
    for script in ${jobscripts}; do
      echo "running '${script}'.."
      ${script} > ${scriptlogprefix}.out
    done
    ;;
  "slurm")
    scriptcats="$(echo ${jobscripts} | grep -o "script[0-9]\+" | sort | uniq )"
    [ ! -z "$scriptcats" ] || { echo "error: no script categories"; exit 1; } >&2
    jobdep=""
    for cat in ${scriptcats}; do
      echo "$cat"
      joblist=""
      for script in $(echo "${jobscripts}" | grep "${cat}_" ); do
        jobname=$( echo $( basename $script ) | grep -o 'script[^/]\+' | sed 's/cript//' )
        echo "${cfg_queuecmd} ${jobdep}
              --job-name=${jobname}
              --output=${tmpprefix}_slurm_%x_%j.out
              --parsable ${script}" | printlog 2
        jobid=$( \
          ${cfg_queuecmd} ${jobdep} \
            --job-name=${jobname} \
            --output=${scriptlogprefix}_slurm_%x_%j.out \
            --parsable ${script}
        )
        joblist="${joblist}:${jobid}"
      done
      jobdep="--dependency=afterok${joblist}"
    done
    ;;
  *)
    echo "error: unknown execmode '${cfg_execmode}'" >&2
    exit 1
    ;;
esac

