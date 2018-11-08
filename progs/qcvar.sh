#!/usr/bin/env bash

# exit on error
trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_freqstd=$( cfgvar_get freqstd )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare    cfg_hweneglogp=$( cfgvar_get hweneglogp )
declare    cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_phenotypes=$( cfgvar_get phenotypes )
declare -r cfg_varmiss=$( cfgvar_get varmiss )

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


# set Hardy-Weinberg test p-value threshold in case of no phenotypic information
[ "${cfg_phenotypes}" != "" -a -s "${cfg_phenotypes}" ] || cfg_hweneglogp=${cfg_hweneglogp_ctrl}

tmp_varmiss=${cfg_varmiss}
n=$( wc -l ${opt_inprefix}.fam | cut -d ' ' -f 1 )
if [ $n -lt ${cfg_minindcount} ] ; then tmp_varmiss=0.1 ; fi
plink --bfile ${opt_inprefix} \
      --not-chr 23,24 \
      --set-hh-missing \
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
sex_hweneglogp=$(( cfg_hweneglogp/2 ))
if [ $sex_hweneglogp -lt 12 ] ; then
  sex_hweneglogp=12
fi
plink --bfile ${opt_inprefix} ${nosex_flag} \
      --chr 23,24 \
      --set-hh-missing \
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

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
mv ${tmpprefix}_out.log ${opt_outprefix}.log

# rm ${tmpprefix}*

