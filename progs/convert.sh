#!/usr/bin/env bash

set -Eeou pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

bcftoolsexec=${BASEDIR}/lib/3rd/bcftools

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_chromosomes=$( cfgvar_get chromosomes )

#-------------------------------------------------------------------------------

# input: final QC plink file set
# output: chromosome-wise eagle VCF input files

printf "\
  * Purge individual and family IDs
  * Convert plink file set into VCF format
" | printlog 1

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi
# purge fam file ids of any unwanted characters
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
               --make-bed --out ${tmpprefix}_idfix \
               2>&1 | printlog 2
else
  printf "====> warning: could not purge IDs due to ensuing conflicts.\n"
  cp ${opt_inprefix}.bed ${tmpprefix}_idfix.bed
  cp ${opt_inprefix}.bim ${tmpprefix}_idfix.bim
  cp ${opt_inprefix}.fam ${tmpprefix}_idfix.fam
fi
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_idfix.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_idfix.fam
declare -ra cfg_chromosomes_arr=( $cfg_chromosomes )
tmp_chromosomes=$( join \
  <( printf "%s\n" ${cfg_chromosomes_arr[@]} | sort -u ) \
  <( cut -f 1 ${tmpprefix}_idfix.bim | sort -u ) \
)
# convert input files to vcf formats
for chr in ${tmp_chromosomes} ; do
  if [ -e ${opt_outprefix}_chr${chr}.bcf ] ; then
    printf "bcf file '%s' found. skipping conversion..\\n" ${opt_outprefix}_chr${chr}.bcf
    continue
  fi
  ${plinkexec} --bfile ${tmpprefix}_idfix \
               --chr ${chr} \
               --recode vcf bgz \
               --out ${tmpprefix}_chr${chr} \
               2>&1 | printlog 2
  ${bcftoolsexec} convert -Ob --threads 3 \
                  ${tmpprefix}_chr${chr}.vcf.gz \
                  > ${tmpprefix}_chr${chr}_out.bcf
  ${bcftoolsexec} annotate --rename-chrs <( echo -e '23 X\n25 X' ) -Ob --threads 3 \
                  ${tmpprefix}_chr${chr}_out.bcf \
                  > ${tmpprefix}_chr${chr}_outx.bcf
  ${bcftoolsexec} index ${tmpprefix}_chr${chr}_outx.bcf
  mv ${tmpprefix}_chr${chr}_outx.bcf     ${opt_outprefix}_chr${chr}.bcf
  mv ${tmpprefix}_chr${chr}_outx.bcf.csi ${opt_outprefix}_chr${chr}.bcf.csi
done

rm -r ${tmpprefix}*

