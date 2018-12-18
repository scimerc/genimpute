#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_freqstd=$( cfgvar_get freqstd )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare    cfg_hweneglogp=$( cfgvar_get hweneglogp )
declare    cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_metrios=$( cfgvar_get metrios )
declare -r cfg_mevars=$( cfgvar_get mevars )
declare -r cfg_phenotypes=$( cfgvar_get phenotypes )
declare -r cfg_varmiss=$( cfgvar_get varmiss )

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping variant QC step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

# input: clean plink set
# output: plink genotype set for variants passing QC criteria:
#         - good coverage
#         - sufficient minor allele frequency
#         - HW equilibrium (possibly different for sex and non-sex chromosomes)
#         - low rate of Mendel errors in trios


declare keepflag=''
# set keep flag if a list of unrelated individuals exists
if [ -f "${opt_outprefixbase}.ids" ] ; then
  keepflag="--keep ${opt_outprefixbase}.ids"
fi

# set Hardy-Weinberg test p-value threshold in case of no phenotypic information
[ "${cfg_phenotypes}" != "" -a -s "${cfg_phenotypes}" ] || cfg_hweneglogp=${cfg_hweneglogp_ctrl}

tmp_varmiss=${cfg_varmiss}
n=$( wc -l ${opt_inprefix}.fam | cut -d ' ' -f 1 )
if [ $n -lt ${cfg_minindcount} ] ; then tmp_varmiss=0.1 ; fi
${plinkexec} --bfile ${opt_inprefix} ${keepflag} \
             --not-chr 23,24 \
             --set-hh-missing \
             --geno ${tmp_varmiss} \
             --maf ${cfg_freqstd} \
             --hwe 1.E-${cfg_hweneglogp} ${cfg_hweflag} \
             --me ${cfg_metrios} ${cfg_mevars} \
             --make-just-bim \
             --out ${tmpprefix}_nonsex \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}
declare nosexflag=''
awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' ${opt_hqprefix}.fam \
  > ${opt_hqprefix}.nosex
if [ -s "${opt_hqprefix}.nosex" ] ; then
  nosexflag="--remove ${opt_hqprefix}.nosex"
fi
sex_hweneglogp=$(( cfg_hweneglogp*2 ))
if [ $sex_hweneglogp -gt 12 ] ; then
  sex_hweneglogp=12
fi
${plinkexec} --bfile ${opt_inprefix} ${keepflag} ${nosexflag} \
             --chr 23,24 \
             --set-hh-missing \
             --geno ${tmp_varmiss} \
             --maf ${cfg_freqstd} \
             --hwe 1.E-${sex_hweneglogp} ${cfg_hweflag} \
             --make-just-bim \
             --out ${tmpprefix}_sex \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}
unset nosexflag
[ -s "${tmpprefix}_nonsex.bim" -o -s "${tmpprefix}_sex.bim" ] || {
  printf "error: no variants left after QC.\n" >&2;
  exit 1;
}
cut -f 2 ${tmpprefix}_*sex.bim | sort -u > ${tmpprefix}.mrk
plink --bfile ${opt_inprefix} \
      --extract ${tmpprefix}.mrk \
      --make-bed \
      --out ${tmpprefix}_out \
      2>&1 >> ${debuglogfn} \
      | tee -a ${debuglogfn}
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
mv ${tmpprefix}_out.log ${opt_outprefix}.log

rm ${tmpprefix}*

