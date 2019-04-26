#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r tmpprefix="${opt_outprefix}_tmp"
declare -r debuglogfn="${tmpprefix}_debug.log"

declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare    cfg_hweneglogp=$( cfgvar_get hweneglogp )
declare    cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_metrios=$( cfgvar_get metrios )
declare -r cfg_mevars=$( cfgvar_get mevars )
declare -r cfg_phenotypes=$( cfgvar_get phenotypes )
declare -r cfg_varmiss=$( cfgvar_get varmiss )

#-------------------------------------------------------------------------------

# input: clean plink set
# output: plink genotype set for variants passing QC criteria:
#         - good coverage
#         - HW equilibrium (possibly different for sex and non-sex chromosomes)
#         - low rate of Mendel errors in trios

printf "\
  * Exclude variants with:
    * low coverage
    * Hardy-Weinberg disequilibrium (separating sex and non-sex chromosomes)
    * high rate of Mendel errors in eventual trios
" | printlog 1

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "'%s' found. skipping variant QC..\n" "${opt_outprefix}.bed"
  exit 0
fi

if ls "${tmpprefix}"* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

printf "removing mendel errors..\n"
${plinkexec} --bfile "${opt_inprefix}" \
             --me ${cfg_metrios} ${cfg_mevars} \
             --set-me-missing \
             --make-bed \
             --out "${tmpprefix}_nome" \
             2> >( tee "${tmpprefix}.err" ) | printlog 2
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi

declare keepflag=''
# set keep flag if a list of unrelated individuals exists
if [ -f "${opt_outprefixbase}.ids" ] ; then
  keepflag="--keep ${opt_outprefixbase}.ids"
fi

# set Hardy-Weinberg test p-value threshold in case of no phenotypic information
[ "${cfg_phenotypes}" != "" -a -s "${cfg_phenotypes}" ] || cfg_hweneglogp=${cfg_hweneglogp_ctrl}

tmp_varmiss=${cfg_varmiss}
n=$( wc -l "${opt_inprefix}.fam" | cut -d ' ' -f 1 )
if [ $n -lt ${cfg_minindcount} ] ; then tmp_varmiss=0.1 ; fi
printf "qc'ing non-sex chromosomes variants..\n"
${plinkexec} --bfile "${tmpprefix}_nome" ${keepflag} \
             --not-chr 23,24 \
             --set-hh-missing \
             --geno ${tmp_varmiss} \
             --hwe 1.E-${cfg_hweneglogp} ${cfg_hweflag} \
             --make-just-bim \
             --out "${tmpprefix}_nonsex" \
             2> >( tee "${tmpprefix}.err" ) | printlog 2
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
declare nosexflag=''
awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' "${opt_hqprefix}.fam" \
  > "${opt_hqprefix}.nosex"
if [ -s "${opt_hqprefix}.nosex" ] ; then
  nosexflag="--remove ${opt_hqprefix}.nosex"
fi
sex_hweneglogp=$(( cfg_hweneglogp*2 ))
if [ $sex_hweneglogp -gt 12 ] ; then
  sex_hweneglogp=12
fi
if [ $( get_xvar_count "${opt_inprefix}.bim" ) -ge ${cfg_minvarcount} ] ; then
  printf "qc'ing sex chromosomes variants..\n"
  ${plinkexec} --bfile "${tmpprefix}_nome" ${keepflag} ${nosexflag} \
               --chr 23,24 \
               --set-hh-missing \
               --geno ${tmp_varmiss} \
               --hwe 1.E-${sex_hweneglogp} ${cfg_hweflag} \
               --make-just-bim \
               --out "${tmpprefix}_sex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 2
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
fi
unset nosexflag
[ -s "${tmpprefix}_nonsex.bim" -o -s "${tmpprefix}_sex.bim" ] || {
  printf "error: no variants left after QC.\n" >&2;
  exit 1;
}
cut -f 2 "${tmpprefix}_"*sex.bim | sort -u > "${tmpprefix}.mrk"
${plinkexec} --bfile "${opt_inprefix}" \
             --extract "${tmpprefix}.mrk" \
             --set-hh-missing \
             --make-bed \
             --out "${tmpprefix}_out" \
             2> >( tee "${tmpprefix}.err" ) | printlog 2
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"

rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out".*

rm "${tmpprefix}"*

