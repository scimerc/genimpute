#!/usr/bin/env bash

set -Eeou pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

execmode="slurm"
# this time is for imputation - expected time for phasing jobs: ~10hrs
slurmcmd="sbatch --account=p33 --time=4-00:00:00 --mem-per-cpu=2048M --cpus-per-task=4" #1536

genmap=${BASEDIR}/lib/data/genetic_map_hg19_withX.txt.gz

inprefix=/cluster/projects/p33/data/genetics/rawdata/genotypes/postQC/NORMENT/norment_batch3_Jan15_qc
outprefix=/cluster/projects/p33/nobackup/tmp/test_batch3
tmpprefix=${outprefix}_tmp
phaserefprefix=/cluster/projects/p33/users/franbe/norment_2018/ega.grch37.chr
imputerefprefix=/cluster/projects/p33/data/genetics/external/HRC/decrypted/ega.grch37.chr
scriptprefix=${tmpprefix}_script

chromosomes=$( seq 13 19 )
igroupsize=650 # num_samples*22 / (400[=max_jobs] - 3[=num_of_chr-wise-jobs]*22[=num_of_chr])
eaglexec=${BASEDIR}/lib/3rd/eagle
minimacexec=${BASEDIR}/lib/3rd/minimac3

getbpstart() {
  local -r chr=$1
  local -r bim=$2
  sort -t $'\t' -k 1,1 ${bim} \
    | join -t $'\t' <( echo $chr ) - \
    | cut -f 4 \
    | sort -n \
    | head -n 1
}

getbpend() {
  local -r chr=$1
  local -r ifn=$2
  sort -t $'\t' -k 1,1 ${bim} \
    | join -t $'\t' <( echo $chr ) - \
    | cut -f 4 \
    | sort -nr \
    | head -n 1
}

if [ ! -s ${tmpprefix}.bcf -o ! -s ${tmpprefix}.bcf.csi ] ; then
  module load plink
  plink --bfile ${inprefix} --recode vcf bgz --out ${tmpprefix}
  bcftools convert -Ob ${tmpprefix}.vcf.gz > ${tmpprefix}_out.bcf
  bcftools index ${tmpprefix}_out.bcf
  mv ${tmpprefix}_out.bcf ${tmpprefix}.bcf
  mv ${tmpprefix}_out.bcf.csi ${tmpprefix}.bcf.csi
fi

#TODO: remove QC'd away samples from bio before splitting
cut -f 1 ${inprefix}.bio | tail -n +2 | split -d -l ${igroupsize} /dev/stdin ${tmpprefix}_sample
tmpsamplelist=$( ls "${tmpprefix}_sample"* | grep 'sample[0-9]\+$' )

for chr in ${chromosomes} ; do
cat > ${scriptprefix}1_phase_chr${chr}.sh << EOI
#!/usr/bin/env bash
set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_chr${chr}_phased.vcf.gz" ]; then
  printf "phased haplotypes present - skipping...\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
  time ${eaglexec} \\
    --chrom ${chr} \\
    --geneticMapFile ${genmap} \\
    --vcfRef ${phaserefprefix}${chr}.haplotypes.bcf \\
    --vcfTarget ${tmpprefix}.bcf \\
    --outPrefix ${tmpprefix}_chr${chr}_phasing \\
    --numThreads \${num_cpus}
mv ${tmpprefix}_chr${chr}_phasing.vcf.gz ${tmpprefix}_chr${chr}_phased.vcf.gz
EOI
done

#TODO? remove me
for chr in ${chromosomes} ; do
cat > ${scriptprefix}2_refconv_chr${chr}.sh << EOI
#!/usr/bin/env bash
set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_refhaps_chr${chr}.m3vcf.gz" ]; then
  printf "m3vcf files present - skipping...\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
    time ${minimacexec} \\ 
      --refHaps ${imputerefprefix}${chr}.haplotypes.vcf.gz \\
      --processReference \\
      --prefix ${tmpprefix}_refhaps_chr${chr}_out \\
      --cpus \${num_cpus} \\
      --rounds 5
    mv ${tmpprefix}_refhaps_chr${chr}_out.m3vcf.gz \\
       ${tmpprefix}_refhaps_chr${chr}.m3vcf.gz \\
EOI
done

for chr in ${chromosomes} ; do
  for samplefile in ${tmpsamplelist} ; do
    sample=${samplefile##*_}
cat > ${scriptprefix}3_impute_chr${chr}_${sample}.sh << EOI
#!/usr/bin/env bash
set -Eeou pipefail
source /cluster/bin/jobsetup
if [ -e "${tmpprefix}_chr${chr}_${sample}_imputed.dose.vcf.gz" ]; then
  printf "m3vcf files present - skipping...\\n"
  exit 0
fi
num_cpus_detected=\$(cat /proc/cpuinfo | grep "model name" | wc -l)
num_cpus=\${OMP_NUM_THREADS:-\${num_cpus_detected}}
    bcftools view -S ${samplefile} -Oz \\
      --force-samples ${tmpprefix}_chr${chr}_phased.vcf.gz \\
      > ${tmpprefix}_chr${chr}_${sample}_phased.vcf.gz
    time ${minimacexec} \\
      --cpus \${num_cpus} \\
      --haps ${tmpprefix}_chr${chr}_${sample}_phased.vcf.gz \\
      --refHaps ${tmpprefix}_refhaps_chr${chr}.m3vcf.gz \\
      --rounds 5 --states 200 \\
      --noPhoneHome \\
      --prefix ${tmpprefix}_chr${chr}_${sample}_imputing \\
      > ${tmpprefix}_chr${chr}_${sample}_imputing.log 2>&1
    bcftools index ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz
    # fix minimac3 output
    bcftools view -h ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz \\
      | awk -v genofilter='##FILTER=<ID=GENOTYPED,Description="Site was genotyped">' \\
        'BEGIN{ flag = 0 } {
          if ( \$0 ~ "^##FILTER" ) {
            if ( flag==0 ) print( genofilter )
            flag = 1
          }
          print
        }' \\
      > ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch
    bcftools reheader -h ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch \\
      ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz \\
      > ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch.dose.vcf.gz
    mv \\
      ${tmpprefix}_chr${chr}_${sample}_imputing_hdrpatch.dose.vcf.gz \\
      ${tmpprefix}_chr${chr}_${sample}_imputed.dose.vcf.gz
    bcftools index -f ${tmpprefix}_chr${chr}_${sample}_imputed.dose.vcf.gz
EOI
  done
cat > ${scriptprefix}4_merge_chr${chr}.sh << EOI
#!/usr/bin/env bash
set -Eeou pipefail
source /cluster/bin/jobsetup
  bcftools merge ${tmpprefix}_chr${chr}_*_imputed.dose.vcf.gz -Oz \\
    > ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz
  mv ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz \\
    ${outprefix}_chr${chr}_imputed.dose.vcf.gz
EOI
done

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

