#!/usr/bin/env bash

set -Eeou pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

execmode="slurm"
slurmcmd="sbatch --account=p33"

genmap=${BASEDIR}/lib/data/genetic_map_hg19_withX.txt.gz

inprefix=/cluster/projects/p33/data/genetics/rawdata/genotypes/postQC/NORMENT/norment_batch3_Jan15_qc
outprefix=/cluster/projects/p33/nobackup/tmp/fk/batch3/batch3
tmpprefix=${outprefix}_tmp
phaserefprefix=/cluster/projects/p33/users/franbe/norment_2018/ega.grch37.chr
imputerefprefix=/cluster/projects/p33/data/genetics/external/HRC/decrypted/ega.grch37.chr
refprefix=/cluster/projects/p33/nobackup/tmp/ref/hrc_refhaps
scriptprefix=${tmpprefix}_script

chromosomes=$( seq 1 22 )
igroupsize=650 # TODO: auto-compute - num_samples*22 / (400[=max_jobs] - 3[=num_of_chr-wise-jobs]*22[=num_of_chr])
eaglexec=${BASEDIR}/lib/3rd/eagle
minimacexec=${BASEDIR}/lib/3rd/Minimac3_march-sb_omp
bcftoolsexec=${BASEDIR}/lib/3rd/bcftools
plinkexec=${BASEDIR}/lib/3rd/plink
timexec='/usr/bin/time -o /dev/stdout -f "ResStats\tRealSec:\t%e\tMaxMemKB:\t%M\tUserSec:\t%U\tSysSec:\t%S"'

#-------------------------------------------------------------------------------

# convert

cat > ${scriptprefix}0_conv.sh << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00

set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e ${tmpprefix}.bcf.csi ] ; then
  printf "input files ready for phasing - skipping...\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
${timexec} ${plinkexec} \\
  --threads \${num_cpus} \\
  --memory 30000 \\
  --bfile ${inprefix} \\
  --recode vcf bgz \\
  --out ${tmpprefix}
${bcftoolsexec} convert -O b ${tmpprefix}.vcf.gz > ${tmpprefix}_out.bcf
${bcftoolsexec} index ${tmpprefix}_out.bcf
mv ${tmpprefix}_out.bcf     ${tmpprefix}.bcf
mv ${tmpprefix}_out.bcf.csi ${tmpprefix}.bcf.csi
EOI

#-------------------------------------------------------------------------------

# phase

#TODO: remove QC'd away samples from bio before splitting
cut -f 1 ${inprefix}.bio | tail -n +2 | split -d -l ${igroupsize} /dev/stdin ${tmpprefix}_sample
tmpsamplelist=$( ls "${tmpprefix}_sample"* | grep 'sample[0-9]\+$' )

for chr in ${chromosomes} ; do
cat > ${scriptprefix}1_phase_chr${chr}.sh << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00

# Eagle's parallelization is very efficient. colossus machines seem to 
# have 20 cpus but 40 hyperthreading cores. Therefore we use 
# num_cpus*2 threads.
# Running time: From 18min to 3h for batch3 (9.9k samples)

set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_chr${chr}_phased.vcf.gz" ]; then
  printf "phased haplotypes present - skipping...\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
${timexec} ${eaglexec} \\
  --chrom ${chr} \\
  --geneticMapFile ${genmap} \\
  --vcfRef ${phaserefprefix}${chr}.haplotypes.bcf \\
  --vcfTarget ${tmpprefix}.bcf \\
  --outPrefix ${tmpprefix}_chr${chr}_phasing \\
  --numThreads \$(( num_cpus * 2))
mv ${tmpprefix}_chr${chr}_phasing.vcf.gz ${tmpprefix}_chr${chr}_phased.vcf.gz
EOI
done

#-------------------------------------------------------------------------------

# impute

for chr in ${chromosomes} ; do
  for samplefile in ${tmpsamplelist} ; do
    sample=${samplefile##*_}
cat > ${scriptprefix}3_impute_chr${chr}_${sample}.sh << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00

# Minimac3-omp cpus=2, mem-per-cpu=16 + lowMem flag works fine with hrc-32k
# (23G mem usage and 300min wall clock time for a chr2 n=650 sample) 
# Number of samples does not seem to increase mem usage.
# Mem usage seems to be determined by ref file.
#
# TODO:try cpu=4 and mem-per-cpu=8 + lowMem flag

set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_chr${chr}_${sample}_imputed.dose.vcf.gz" ]; then
  printf "m3vcf files present - skipping...\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
${bcftoolsexec} view -S ${samplefile} -Oz \\
  --force-samples ${tmpprefix}_chr${chr}_phased.vcf.gz \\
  > ${tmpprefix}_chr${chr}_${sample}_phased.vcf.gz
${timexec} ${minimacexec} \\
  --cpus \${num_cpus} \\
  --lowMemory \\
  --haps ${tmpprefix}_chr${chr}_${sample}_phased.vcf.gz \\
  --refHaps ${refprefix}_chr${chr}.m3vcf.gz \\
  --rounds 5 --states 200 \\
  --noPhoneHome \\
  --prefix ${tmpprefix}_chr${chr}_${sample}_imputing \\
  > ${tmpprefix}_chr${chr}_${sample}_imputing.log 2>&1
${bcftoolsexec} index ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz
# fix minimac3 output
${bcftoolsexec} view -h ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz \\
  | awk -v genofilter='##FILTER=<ID=GENOTYPED,Description="Site was genotyped">' \\
    'BEGIN{ flag = 0 } {
      if ( \$0 ~ "^##FILTER" ) {
        if ( flag==0 ) print( genofilter )
        flag = 1
      }
      print
    }' \\
  > ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch
${bcftoolsexec} reheader -h ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch \\
  ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz \\
  > ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch.dose.vcf.gz
mv \\
  ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch.dose.vcf.gz \\
  ${tmpprefix}_chr${chr}_${sample}_imputed.dose.vcf.gz
${bcftoolsexec} index -f ${tmpprefix}_chr${chr}_${sample}_imputed.dose.vcf.gz
EOI
  done
done

#-------------------------------------------------------------------------------

# merge

for chr in ${chromosomes} ; do
  cat > ${scriptprefix}4_merge_chr${chr}.sh << EOI
#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00

set -Eeou pipefail
source /cluster/bin/jobsetup
${bcftoolsexec} merge ${tmpprefix}_chr${chr}_*_imputed.dose.vcf.gz -Oz \\
  > ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz
mv ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz \\
   ${outprefix}_chr${chr}_imputed.dose.vcf.gz
EOI
done

#-------------------------------------------------------------------------------

# exec

echo 'submitting scripts..'

jobscripts=$(ls ${scriptprefix}*)

case ${execmode} in
  "local")
    for script in ${jobscripts}; do
      echo "running" ${script} "..."
      ${script}
    done
    ;;
  "slurm")
    scriptcats="$(echo ${jobscripts} | grep -o "script[0-9]\+" | sort | uniq )"
    [ ! -z "$scriptcats" ] || { printf "error: no script categories"; exit 1; } >&2
    jobdep=""
    for cat in ${scriptcats}; do
      echo "$cat"
      joblist=""
      for script in $(echo "${jobscripts}" | grep "_${cat}_" ); do
        jobname=$( echo $script | grep -o 'script.\+' | sed 's/cript//' )
        jobid=$( \
          ${slurmcmd} ${jobdep} \
            --job-name=${jobname} \
            --output=${tmpprefix}_slurm_%x_%j.out \
            --parsable ${script}
        )
        joblist="${joblist}:${jobid}"
      done
      jobdep="--dependency=afterok${joblist}"
    done
    ;;
  *)
    printf "error unknown execmode '%s'\n" "${execmode}" >&2
    exit 1
    ;;
esac

