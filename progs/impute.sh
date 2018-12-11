#!/usr/bin/env bash

set -Eeou pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

genmap=${BASEDIR}/lib/data/genetic_map_hg19_withX.txt.gz

inprefix=/cluster/projects/p33/nobackup/tmp/testZ
outprefix=/cluster/projects/p33/nobackup/tmp/testZQ
tmpprefix=${outprefix}_tmp
phaserefprefix=/cluster/projects/p33/users/franbe/norment_2018/test/chrfiles/chr
imputerefprefix=/cluster/projects/p33/users/franbe/norment_2018/test/chrfiles/chr
scriptprefix=${tmpprefix}_script

igroupsize=100
eaglexec=eagle
minimacexec=minimac

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
  plink --bfile ${inprefix} --recode vcf bgz --out ${tmpprefix}
  bcftools convert -Ob ${tmpprefix}.vcf.gz > ${tmpprefix}_out.bcf
  mv ${tmpprefix}_out.bcf ${tmpprefix}.bcf
  bcftools index ${tmpprefix}.bcf
fi

cut -f 1 ${inprefix}.bio | tail -n +2 | split -d -l ${igroupsize} /dev/stdin ${tmpprefix}_sample
tmpsamplelist=$( ls "${tmpprefix}_sample"* | grep 'sample[0-9]\+$' )

for chr in $( seq 21 22 ) ; do
cat > ${scriptprefix}_phase_chr${chr}.sh << EOI
set -Eeou pipefail
trap 'printf "===> error in %s line %s\n" \$(basename \$0) \${LINENO}; exit;' ERR
num_cpus=$(cat /proc/cpuinfo | grep "model name" | wc -l)
  time ${eaglexec} \\
    --chrom ${chr} \\
    --geneticMapFile ${genmap} \\
    --vcfRef ${phaserefprefix}${chr}ref.vcf.gz \\
    --vcfTarget ${tmpprefix}.bcf \\
    --outPrefix ${tmpprefix}_chr${chr}_phasing \\
    --numThreads \${num_cpus}
mv ${tmpprefix}_chr${chr}_phasing.vcf.gz ${tmpprefix}_chr${chr}_phased.vcf.gz
EOI
done

for chr in $( seq 21 22 ) ; do
  for samplefile in ${tmpsamplelist} ; do
    sample=${samplefile##*_}
cat > ${scriptprefix}_impute_chr${chr}_${sample}.sh << EOI
set -Eeou pipefail
trap 'printf "===> error in %s line %s\n" \$(basename \$0) \${LINENO}; exit;' ERR
num_cpus=$(cat /proc/cpuinfo | grep "model name" | wc -l)
    bcftools view -S ${samplefile} -Oz \\
      --force-samples ${tmpprefix}_chr${chr}_phased.vcf.gz \\
      > ${tmpprefix}_chr${chr}_${sample}_phased.vcf.gz
    time ${minimacexec} \\
      --chr ${chr} \\
      --cpus \${num_cpus} \\
      --haps ${tmpprefix}_chr${chr}_${sample}_phased.vcf.gz \\
      --refHaps ${imputerefprefix}${chr}ref.vcf.gz \\
      --rounds 5 --states 200 \\
      --noPhoneHome --doseOutput --hapOutput \\
      --prefix ${tmpprefix}_chr${chr}_${sample}_imputing \\
      > ${tmpprefix}_chr${chr}_${sample}_imputing.log 2>&1
    bcftools index ${tmpprefix}_chr${chr}_${sample}_imputing.dose.vcf.gz
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
cat > ${scriptprefix}_merge_chr${chr}.sh << EOI
  bcftools merge ${tmpprefix}_chr${chr}_*_imputed.dose.vcf.gz -Oz \\
    > ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz
  mv ${tmpprefix}_chr${chr}_imputed.dose.vcf.gz \\
    ${outprefix}_chr${chr}_imputed.dose.vcf.gz
EOI
done

