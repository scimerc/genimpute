#!/usr/bin/env bash

set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

bcftoolsexec="${BASEDIR}/lib/3rd/bcftools"
eaglexec="${BASEDIR}/lib/3rd/eagle"
minimac_version=3

declare -r cfg_chromosomes=$( cfgvar_get chromosomes )
declare -r cfg_execommand=$( cfgvar_get execommand )
declare -r cfg_genomemap=$( cfgvar_get genomemap )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare -r cfg_hweneglogp=$( cfgvar_get hweneglogp )
declare -r cfg_igroupsize=$( cfgvar_get igroupsize )
declare -r cfg_impmaf=$( cfgvar_get impmaf )
declare -r cfg_impmingp=$( cfgvar_get impmingp )
declare -r cfg_imprsqbase=$( cfgvar_get imprsqbase )
declare -r cfg_imprsqhigh=$( cfgvar_get imprsqhigh )
declare -r cfg_impvarmiss=$( cfgvar_get impvarmiss )
declare -r cfg_metrios=$( cfgvar_get metrios )
declare -r cfg_mevars=$( cfgvar_get mevars )
declare -r genmapfile="${BASEDIR}/lib/data/${cfg_genomemap}"
declare -r scriptlogprefix="${opt_outprefixbase}/.i/.s/logs/script"
declare -r scriptprefix="${opt_outprefixbase}/.i/.s/scripts/script"
declare -r sampleprefix="${opt_outprefixbase}/.i/.s/samples/sample"
declare -r tmpprefix="${opt_outprefix}_tmp"
declare -a jobscripts

mkdir -p "$( dirname "${sampleprefix}" )"
mkdir -p "$( dirname "${scriptprefix}" )"
mkdir -p "$( dirname "${scriptlogprefix}" )"

#-------------------------------------------------------------------------------

# extract list of individuals to establish unique order
for bfile in "${opt_inprefix}_chr"*.bcf ; do
  "${bcftoolsexec}" query -l "${bfile}"
done | sort -u | awk '{ print( $1, $1 ) }' > "${tmpprefix}_ordered.fam"

# split individuals in groups
#TODO: enable sex-specific subdivisions
printf "> splitting individuals into groups..\n"
for bfile in "${opt_inprefix}_chr"*.bcf ; do
  "${bcftoolsexec}" query -l "${bfile}"
done | sort -u | split -d -l ${cfg_igroupsize} /dev/stdin "${sampleprefix}"
groupsize=${cfg_igroupsize}
tmplist=$( ls "${sampleprefix}"* )
Nind=$( cat "${sampleprefix}"* | wc -l )
if [ ${Nind} -lt ${groupsize} ] ; then
  groupsize=${Nind}
fi

#-------------------------------------------------------------------------------

echo -e "==== Eagle phasing ====\n" | printlog 1

# write phase scripts
printf "> writing phasing scripts..\n"
for chr in ${cfg_chromosomes} ; do
  chrtag=${chr}
  if [ ${chr} -eq 23 -o ${chr} -eq 25 ] ; then
    chrtag=X
  fi
  Nvar=$( "${bcftoolsexec}" query -f '%ID\n' "${opt_inprefix}_chr${chr}.bcf" | wc -l )
#TODO: implement different running time for reference-less phasing (stage ii)
  runtimehrs=$(( 1 + Nind / 10000 + Nvar / 2000 ))
  phasescriptfn="${scriptprefix}1_phase_chr${chr}.sh"
  cat > "${phasescriptfn}" << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=512MB
#SBATCH --time=${runtimehrs}:00:00

# Eagle's parallelization is very efficient. The original colossus machines had
# 20 cpus but 40 hyperthreading cores. Therefore, we used num_cpus*2 threads.
# Running time: 10c-2G-20t From 18min to 3h for batch3 (~9k individuals)

set -Eeou pipefail
[ -s /cluster/bin/jobsetup ] && source /cluster/bin/jobsetup
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
${timexec} "${eaglexec}" \\
  --chrom ${chrtag} \\
  --geneticMapFile "${genmapfile}" \\
  --vcfRef "${cfg_refprefix}.chr${chr}.haplotypes.bcf" \\
  --vcfTarget "${opt_inprefix}_chr${chr}.bcf" \\
  --outPrefix "${tmpprefix}_chr${chr}_phasing" \\
  --numThreads \${num_cpus}
mv "${tmpprefix}_chr${chr}_phasing.vcf.gz" "${tmpprefix}_chr${chr}_phased.vcf.gz"
EOI
  chmod u+x ${phasescriptfn}
  if [ -e "${tmpprefix}_chr${chr}_phased.vcf.gz" ]; then
    printf "> phased haplotypes present. nothing to do.\n"
    continue
  fi
  if [ ! -e "${cfg_refprefix}.chr${chr}.haplotypes.bcf" ] ; then
    printf "> reference for chromosome ${chr} not found. skipping phasing..\n"
    continue
  fi
  printf "> adding ${phasescriptfn} to jobs stack..\n"
  jobscripts+=( "${phasescriptfn}" )
done

#-------------------------------------------------------------------------------

echo "==== MaCH${minimac_version} imputation ====" | printlog 1

# write impute scripts
printf "> writing imputation scripts..\n"
SBATCH_CONF_MM3="
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=$(( 12 * ( groupsize / 500 + 1 ) )):00:00

# results for batch3; n=9200 individuals total, 14 samples with n=650 each; hrc-32k reference
# Minimac3-omp
#   Number of individuals does not seem to increase mem usage.
#   Mem usage seems to be determined by ref file.
#   3) sbatch: cpus-per-task=4; mem-per-cpu=14 mm3: lowMem; cpus=8  -> 10h; 54jobs; chr2 max.mem=52/56G
#   2) sbatch: cpus-per-task=4; mem-per-cpu=8  mm3: lowMem; cpus=4  -> 14h; 46jobs; chr2 max.mem=31/32G(!)
#   1) sbatch: cpus-per-task=2; mem-per-cpu=16 mm3: lowMem; cpus=2  -> 21h; 46jobs; chr2 max.mem=23/32G
"
SBATCH_CONF_MM4="
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3GB
#SBATCH --time=$(( 4 * ( groupsize / 500 + 1 ) )):00:00

# results for batch3; n=9200 individuals total, 14 samples with n=650 each; hrc-32k reference
# minimac4
#   5) sbatch: cpus-per-task=4 mem-per-cpu=4G  mm4: cpus=8  -> 1.5h; 178jobs; chr6 max.mem=7/16G
#   4) sbatch: cpus-per-task=2 mem-per-cpu=8G  mm4: cpus=4  -> 3.0h; 126jobs; chr6 max.mem=6/16G
"

for chr in ${cfg_chromosomes} ; do
  for samplefile in ${tmplist} ; do
    sampletag="${samplefile##*\/}"
    sampletag="${sampletag//*sample/s}"
    sampletag="${sampletag//*s/sample}"
    imputescriptfn="${scriptprefix}2_impute_chr${chr}_${sampletag}.sh"
    case "${minimac_version}" in
      "3")
        minimacexec="${BASEDIR}/lib/3rd/minimac3_march-sb_omp --lowMemory"
        sbatch_conf=${SBATCH_CONF_MM3}
        ;;
      "4")
        minimacexec="${BASEDIR}/lib/3rd/minimac4"
        sbatch_conf=${SBATCH_CONF_MM4}
        ;;
      *)
        printf "> error unknown minimac_version '%s'\n" "${minimac_version}" >&2
        exit 1
        ;;
    esac
    cat > "${imputescriptfn}" << EOI
#!/usr/bin/env bash
${sbatch_conf}

set -Eeou pipefail
[ -s /cluster/bin/jobsetup ] && source /cluster/bin/jobsetup
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
"${bcftoolsexec}" view -S "${samplefile}" -Oz --threads 3 \\
  --force-samples "${tmpprefix}_chr${chr}_phased.vcf.gz" \\
  > "${tmpprefix}_chr${chr}_${sampletag}_phased_draft.vcf.gz"
mv "${tmpprefix}_chr${chr}_${sampletag}_phased_draft.vcf.gz" \\
  "${tmpprefix}_chr${chr}_${sampletag}_phased.vcf.gz"
${timexec} ${minimacexec} \\
  --cpus \${num_cpus} \\
  --format GT,GP,DS \\
  --haps "${tmpprefix}_chr${chr}_${sampletag}_phased.vcf.gz" \\
  --refHaps "${cfg_refprefix}.chr${chr}.m3vcf.gz" \\
  --noPhoneHome \\
  --prefix "${tmpprefix}_chr${chr}_${sampletag}_imputing" \\
  > "${tmpprefix}_chr${chr}_${sampletag}_imputing.log" 2>&1
"${bcftoolsexec}" index -f "${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.vcf.gz"
# fix minimac3 output
"${bcftoolsexec}" view -h "${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.vcf.gz" \\
  | awk -v genofilter='##FILTER=<ID=GENOTYPED,Description="Site was genotyped">' \\
    'BEGIN{ flag = 0 } {
      if ( \$0 ~ "^##FILTER" ) {
        if ( flag==0 ) print( genofilter )
        flag = 1
      }
      print
    }' \\
  > "${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch"
"${bcftoolsexec}" reheader -h "${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch" \\
  "${tmpprefix}_chr${chr}_${sampletag}_imputing.dose.vcf.gz" \\
  > "${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch.dose.vcf.gz"
"${bcftoolsexec}" annotate --set-id %CHROM:%POS:%REF:%ALT -Oz --threads 3 \\
  "${tmpprefix}_chr${chr}_${sampletag}_imputing_hdrpatch.dose.vcf.gz" \\
  > "${tmpprefix}_chr${chr}_${sampletag}_imputing_std.dose.vcf.gz"
"${bcftoolsexec}" index -f "${tmpprefix}_chr${chr}_${sampletag}_imputing_std.dose.vcf.gz"
rename imputing_std.dose.vcf.gz imputed.dose.vcf.gz \\
  "${tmpprefix}_chr${chr}_${sampletag}_imputing_std.dose.vcf.gz"*
EOI
    chmod u+x "${imputescriptfn}"
    if [ -e "${tmpprefix}_chr${chr}_${sampletag}_imputed.dose.vcf.gz" ]; then
      printf "> final vcf files present. nothing to do.\n"
      continue
    fi
    if [ ! -e "${cfg_refprefix}.chr${chr}.m3vcf.gz" ] ; then
      printf "> m3vcf reference for chromosome ${chr} not found. skipping imputation..\n"
      continue
    fi
    printf "> adding ${imputescriptfn} to jobs stack..\n"
    jobscripts+=( "${imputescriptfn}" )
  done
done

#-------------------------------------------------------------------------------

# write merge scripts
printf "> writing merge scripts..\n"
for chr in ${cfg_chromosomes} ; do
  scriptfn="${scriptprefix}3_merge_chr${chr}.sh"
  cat > "${scriptfn}" << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=$(( 6*( Nind / 10000 + 1 ) )):00:00

set -Eeou pipefail
[ -s /cluster/bin/jobsetup ] && source /cluster/bin/jobsetup
Ns=\$( ls "${tmpprefix}_chr${chr}"_*_imputed.dose.vcf.gz | wc -l ) 
if [ \${Ns} -eq 1 ] ; then
  cp "${tmpprefix}_chr${chr}"_*_imputed.dose.vcf.gz "${tmpprefix}_chr${chr}_imputed.dose.vcf.gz"
elif [ \${Ns} -gt 1 ] ; then 
  "${bcftoolsexec}" merge "${tmpprefix}_chr${chr}"_*_imputed.dose.vcf.gz --threads 3 -Oz \\
    > "${tmpprefix}_chr${chr}_imputed.dose.vcf.gz"
fi
"${bcftoolsexec}" filter -e "R2 < ${cfg_imprsqbase}" \\
  "${tmpprefix}_chr${chr}_imputed.dose.vcf.gz" --threads 3 -Oz \\
  > "${tmpprefix}_chr${chr}_qc0_imputed.dose.vcf.gz"
"${bcftoolsexec}" index "${tmpprefix}_chr${chr}_qc0_imputed.dose.vcf.gz"
rename "${tmpprefix}_chr${chr}_qc0_imputed.dose.vcf.gz" \\
       "${opt_outprefixbase}/bcf/qc0/chr${chr}.vcf.gz" \\
       "${tmpprefix}_chr${chr}_qc0_imputed.dose.vcf.gz"*
EOI
  chmod u+x "${scriptfn}"
  if [ -e "${opt_outprefixbase}/bcf/qc0/chr${chr}.vcf.gz" ]; then
    printf "> merged vcf files present. nothing to do.\n"
    continue
  fi
  printf "> adding ${scriptfn} to jobs stack..\n"
  jobscripts+=( "${scriptfn}" )
done

#-------------------------------------------------------------------------------

# write post-imputation qc scripts
printf "> writing post-imputation qc scripts..\n"
for chr in ${cfg_chromosomes} ; do
  scriptfn="${scriptprefix}4_pimp_chr${chr}.sh"
  cat > "${scriptfn}" << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=$(( 6*( Nind / 10000 + 1 ) )):00:00

set -Eeou pipefail
[ -s /cluster/bin/jobsetup ] && source /cluster/bin/jobsetup
"${bcftoolsexec}" filter -e "R2 < ${cfg_imprsqhigh} || MAF < ${cfg_impmaf}" \\
  "${opt_outprefixbase}/bcf/qc0/chr${chr}.vcf.gz" --threads 3 -Oz \\
  > "${tmpprefix}_chr${chr}_qc_imputed.dose.vcf.gz"
"${bcftoolsexec}" index "${tmpprefix}_chr${chr}_qc_imputed.dose.vcf.gz"
rename "${tmpprefix}_chr${chr}_qc_imputed.dose.vcf.gz" \\
       "${opt_outprefixbase}/bcf/qc1/chr${chr}.vcf.gz" \\
       "${tmpprefix}_chr${chr}_qc_imputed.dose.vcf.gz"*
EOI
  chmod u+x "${scriptfn}"
  if [ -e "${opt_outprefixbase}/bcf/qc1/chr${chr}.vcf.gz" ]; then
    printf "> qc'd vcf files present. nothing to do.\n"
    continue
  fi
  printf "> adding ${scriptfn} to jobs stack..\n"
  jobscripts+=( "${scriptfn}" )
done

#-------------------------------------------------------------------------------

# write recode scripts
printf "> writing recode scripts..\n"
for chr in ${cfg_chromosomes} ; do
  scriptfn="${scriptprefix}5_recode_chr${chr}.sh"
  cat > "${scriptfn}" << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=${plink_num_cpus}
#SBATCH --mem-per-cpu=${plink_mem_per_cpu}MB
#SBATCH --time=$(( 6*( Nind / 10000 + 1 ) )):00:00

set -Eeou pipefail
[ -s /cluster/bin/jobsetup ] && source /cluster/bin/jobsetup
${plinkexec} --allow-extra-chr \\
  --vcf "${opt_outprefixbase}/bcf/qc1/chr${chr}.vcf.gz" \\
  --vcf-min-gp ${cfg_impmingp} \\
  --double-id \\
  --make-bed \\
  --out "${tmpprefix}_chr${chr}_plink"
declare nosexflag=''
awk '{ OFS="\t"; if ( NR > 1 && \$3 == 0 ) print( \$1, \$2 ); }' "${opt_inprefix}_updatesex.txt" \
  > "${tmpprefix}.nosex"
if [ -s "${tmpprefix}.nosex" -a ${chr} -ge 23 -a ${chr} -lt 25 ] ; then
  nosexflag="--remove ${tmpprefix}.nosex"
fi
${plinkexec} --allow-extra-chr \\
             --bfile "${tmpprefix}_chr${chr}_plink" \\
             --update-sex "${opt_inprefix}_updatesex.txt" \\
             --make-bed \\
             --out "${tmpprefix}_chr${chr}_sup"
${plinkexec} --allow-extra-chr \\
             --bfile "${tmpprefix}_chr${chr}_sup" \\
             --update-ids "${opt_inprefix}_updateids.txt" \\
             --make-bed \\
             --out "${tmpprefix}_chr${chr}_kingids"
${plinkexec} --allow-extra-chr \\
             --bfile "${tmpprefix}_chr${chr}_kingids" \\
             --update-parents "${opt_inprefix}_updateparents.txt" \\
             --make-bed \\
             --out "${tmpprefix}_chr${chr}_kingpeds"
${plinkexec} --allow-extra-chr \\
             --bfile "${tmpprefix}_chr${chr}_kingpeds" \\
             --me ${cfg_metrios} ${cfg_mevars} \\
             --set-me-missing \\
             --make-bed \\
             --out "${tmpprefix}_chr${chr}_nome"
${plinkexec} --allow-extra-chr \${nosexflag} \\
             --bfile "${tmpprefix}_chr${chr}_nome" \\
             --geno ${cfg_impvarmiss} \\
             --hwe 1.E-${cfg_hweneglogp} ${cfg_hweflag} \\
             --make-bed \\
             --out "${tmpprefix}_chr${chr}_out"
unset nosexflag
${plinkexec} --allow-extra-chr \\
  --bfile "${tmpprefix}_chr${chr}_out" \\
  --indiv-sort f "${tmpprefix}_ordered.fam" \\
  --make-bed \\
  --out "${tmpprefix}_chr${chr}_plink_reordered"
rename "${tmpprefix}_chr${chr}_plink_reordered" \\
       "${opt_outprefixbase}/bed/chr${chr}" \\
       "${tmpprefix}_chr${chr}_plink_reordered"*
EOI
  chmod u+x "${scriptfn}"
  if [ -s "${opt_outprefixbase}/bed/chr${chr}.bed" ]; then
    printf "> plink-recoded files present. nothing to do.\n"
    continue
  fi
  printf "> adding ${scriptfn} to jobs stack..\n"
  jobscripts+=( "${scriptfn}" )
done

#-------------------------------------------------------------------------------

# write plink merge scripts
printf "> writing plink merge script..\n"
scriptfn="${scriptprefix}6_chrcat.sh"
cat > "${scriptfn}" << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=${plink_num_cpus}
#SBATCH --mem-per-cpu=${plink_mem_per_cpu}MB
#SBATCH --time=$(( 6*( Nind / 10000 + 1 ) )):00:00

set -Eeou pipefail
[ -s /cluster/bin/jobsetup ] && source /cluster/bin/jobsetup
ls -1 "${opt_outprefixbase}/bed/chr"*.bed | sed 's/\.bed$//g;' > ${tmpprefix}.list
${plinkexec} --allow-extra-chr \\
  --merge-list "${tmpprefix}.list" \\
  --out "${tmpprefix}_all"
rename "${tmpprefix}_all" \\
       "${opt_outprefixbase}/bed/all" \\
       "${tmpprefix}_all"*
EOI
chmod u+x "${scriptfn}"
if [ -s "${opt_outprefixbase}/bed/all.bed" ]; then
  printf "> plink-merged file present. nothing to do.\n"
  continue
fi
printf "> adding ${scriptfn} to jobs stack..\n"
jobscripts+=( "${scriptfn}" )

#-------------------------------------------------------------------------------

# exec

if [ ${opt_dryimpute} -eq 1 ] ; then
  exit 0
fi

echo '> processing scripts..'

case ${cfg_execommand} in
  *bash*)
    > "${scriptlogprefix}.out"
    for script in ${jobscripts[@]}; do
      echo "> running '${script}'.."
      "${script}" >> "${scriptlogprefix}.out"
    done
    ;;
  *sbatch*)
    scriptcats="$(printf "%s\n" ${jobscripts[@]} | grep -o "script[0-9]\+" | sort | uniq )"
    [ ! -z "$scriptcats" ] || { echo "error: no script categories"; exit 1; } >&2
    ijobdep=""
    for cat in ${scriptcats}; do
      ijoblist=""
      echo "> script category '${cat}':" | printlog 3
      for script in $(printf "%s\n" ${jobscripts[@]} | grep "${cat}_" ); do
        ijobname=$( echo $( basename "${script}" ) | grep -o 'script[^/]\+' | sed 's/cript//' )
        echo " > submitting '${cfg_execommand} ${ijobdep}
              --job-name=${ijobname}
              --output=${scriptlogprefix}_slurm_%x_%j.out
              --parsable ${script}'.. " | printlog 3
        ijobid=$( \
          ${cfg_execommand} ${ijobdep} \
          --job-name=${ijobname} \
          --output="${scriptlogprefix}_slurm_%x_%j.out" \
          --parsable "${script}"
        )
        ijoblist="${ijoblist}:${ijobid}"
        echo ' > done.' | printlog 3
      done
      ijobdep="--dependency=afterok${ijoblist}"
      echo '> done.' | printlog 3
    done
    ;;
  *)
    echo "> error: unknown execution command '${cfg_execommand}'" >&2
    exit 1
    ;;
esac

