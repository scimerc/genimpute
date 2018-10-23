#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log
declare -r cfg_freqstd=0.01
declare -r cfg_hweflag='midp include-nonctrl'
declare -r cfg_hweneglogp=12
declare -r cfg_minindcount=100
declare -r cfg_varmiss=0.05

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping variant QC step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
  exit 1
fi

# input: clean plink set
# output: plink genotype set for variants passing QC criteria:
#         - good coverage
#         - sufficient minor allele frequency
#         - HW equilibrium (possibly different for sex and non-sex chromosomes)

tmp_varmiss=${cfg_varmiss}
n=$( wc -l ${opt_inprefix}.fam | cut -d ' ' -f 1 )
if [ $n -lt ${cfg_minindcount} ] ; then tmp_varmiss=0.1 ; fi
sex_hweneglogp=$(( cfg_hweneglogp*2 ))
if [ $sex_hweneglogp -gt 12 ] ; then
  sex_hweneglogp=12
fi
plink --bfile ${opt_inprefix} \
      --not-chr 23,24 \
      --geno ${tmp_varmiss} \
      --maf ${cfg_freqstd} \
      --hwe 1.E-${cfg_hweneglogp} ${cfg_hweflag} \
      --make-just-bim \
      --out ${tmpprefix}_nonsex \
      >> ${debuglogfn}
nosex_flag=''
awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' ${opt_hqprefix}.fam \
  > ${opt_hqprefix}.nosex
if [ -s "${opt_hqprefix}.nosex" ] ; then
  nosex_flag="--remove ${opt_hqprefix}.nosex"
fi
plink --bfile ${opt_inprefix} ${nosex_flag} \
      --chr 23,24 \
      --geno ${tmp_varmiss} \
      --maf ${cfg_freqstd} \
      --hwe 1.E-${sex_hweneglogp} ${cfg_hweflag} \
      --make-just-bim \
      --out ${tmpprefix}_sex \
      >> ${debuglogfn}
[ -s "${tmpprefix}_nonsex.bim" -o -s "${tmpprefix}_sex.bim" ] || {
  printf "error: no variants left after QC.\n" >&2;
  exit 1;
}
cut -f 2 ${tmpprefix}_*sex.bim | sort -u > ${tmpprefix}.mrk
plink --bfile ${opt_inprefix} \
      --extract ${tmpprefix}.mrk \
      --make-bed \
      --out ${tmpprefix}_out \
      >> ${debuglogfn}
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam

mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.bed ${opt_outprefix}.bed

rm ${tmpprefix}*

