#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r tmpprefix=${opt_hqprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_varmiss=$( cfgvar_get varmiss )
declare -r cfg_freqhq=$( cfgvar_get freqhq )
declare -r cfg_genomeblacklist=$( cfgvar_get genomeblacklist )
declare -r cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_uid=$( cfgvar_get uid )

if [ -f "${opt_hqprefix}.bed" -a -f "${opt_hqprefix}.bim" -a -f "${opt_hqprefix}.fam" ] ; then
  printf "skipping hq step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1 
fi

# input: merged plink set
# output: hq plink set (with imputed sex, if possible)
# 1) get sex hq-variants from input file
# 2) get non-sex hq-variants from input file
# 3) extract all hq-variants from input file and make hq plink set
# 4) LD-prune hq variants from 3)
# 5) impute sex once with all standard hq-variants from 4
# 6) if sex could be imputed for enough individuals, then
#    impute it once again after HWE tests


declare keepflag=''
# set keep flag if a list of unrelated individuals exists
if [ -f "${opt_outprefixbase}.ids" ] ; then
  keepflag="--keep ${opt_outprefixbase}.ids"
fi
declare extractflag=''
# set extract flag if a previous list of variants exists
if [ -f "${opt_outprefixbase}.mrk" ] ; then
  extractflag="--extract ${opt_outprefixbase}.mrk"
fi

declare -r regionblacklist=${BASEDIR}/lib/data/${cfg_genomeblacklist}
# check if exclude file exists and is not empty
[ -s "${regionblacklist}" ] || {
  printf "error: file '%s' empty or not found.\n" ${regionblacklist} >&2;
  exit 1;
}
declare -r regexcludeflag="--exclude range ${regionblacklist}"

# get non sex hq-variants from input file
${plinkexec} --bfile ${opt_inprefix} ${keepflag} \
             --not-chr 23,24 ${regexcludeflag} ${extractflag} \
             --geno ${cfg_varmiss} \
             --maf ${cfg_freqhq} \
             --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag} \
             --make-just-bim \
             --out ${tmpprefix}_nonsex \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}

if [ $( get_xvar_count ${opt_inprefix}.bim ) -ge ${cfg_minvarcount} ] ; then
  # get sex hq-variants from input file
  ${plinkexec} --bfile ${opt_inprefix} ${keepflag} \
               --chr 23,24 ${regexcludeflag} ${extractflag} \
               --geno ${cfg_varmiss} \
               --maf ${cfg_freqhq} \
               --make-just-bim \
               --out ${tmpprefix}_sex \
               2>&1 >> ${debuglogfn} \
               | tee -a ${debuglogfn}
fi

# check if we have anything of high quality
[ -s "${tmpprefix}_nonsex.bim" -o -s "${tmpprefix}_sex.bim" ] || {
  printf "error: no variants left in high quality set." >&2;
  exit 1;
}

# extract all hq variants from input file and make hq plink set
cut -f 2 ${tmpprefix}_*sex.bim | sort -u > ${tmpprefix}_hq.mrk
${plinkexec} --bfile ${opt_inprefix} \
             --extract ${tmpprefix}_hq.mrk \
             --make-bed \
             --out ${tmpprefix}_hq \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}

# LD-prune hq variants
${plinkexec} --bfile ${tmpprefix}_hq ${keepflag} \
             --indep-pairphase 500 5 0.2 \
             --out ${tmpprefix}_hq_LD \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}

# extract LD-pruned hq variants from hq plink set
${plinkexec} --bfile ${tmpprefix}_hq \
             --extract ${tmpprefix}_hq_LD.prune.in \
             --make-bed \
             --out ${tmpprefix}_hq_LDpruned \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}

# if there are enough X chromosome variants impute sex based on them
if [ $( get_xvar_count ${tmpprefix}_hq_LDpruned.bim ) -ge $cfg_minvarcount ] ; then
  # impute sex once with all standard high quality variants
  ${plinkexec} --bfile ${tmpprefix}_hq_LDpruned \
               --impute-sex \
               --make-bed \
               --out ${tmpprefix}_hq_LDpruned_isex \
               2>&1 >> ${debuglogfn} \
               | tee -a ${debuglogfn}

  rename hq_LDpruned_isex out ${tmpprefix}_hq_LDpruned_isex.*

  declare -r xindcount=$( awk '$5 == 1 || $5 == 2' ${tmpprefix}_out.fam | wc -l )
  # if sex could be imputed for enough individuals impute it once again after HWE tests
  if [ ${xindcount} -gt ${cfg_minindcount} ] ; then
    ${plinkexec} --bfile ${tmpprefix}_out \
                 --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag} \
                 --make-just-bim \
                 --out ${tmpprefix}_sexhwe \
                 2>&1 >> ${debuglogfn} \
                 | tee -a ${debuglogfn}

    # if there are enough X chromosome variants after HWE re-impute sex based on them
    if [ $( get_xvar_count ${tmpprefix}_sexhwe.bim ) -ge ${cfg_minvarcount} ] ; then
      ${plinkexec} --bfile ${tmpprefix}_hq_LDpruned \
                   --extract <( cut -f 2 ${tmpprefix}_sexhwe.bim ) \
                   --impute-sex \
                   --make-bed \
                   --out ${tmpprefix}_hq_LDpruned_isex_new \
                   2>&1 >> ${debuglogfn} \
                   | tee -a ${debuglogfn}

      # replace the original sex imputation files
      rename hq_LDpruned_isex_new out ${tmpprefix}_hq_LDpruned_isex_new.*

    fi
  fi

else
  
  # if sex could not be imputed use LD-pruned set
  rename hq_LDpruned out ${tmpprefix}_hq_LDpruned.*

fi

rename ${tmpprefix}_out ${opt_hqprefix} ${tmpprefix}_out.*

rm -f ${tmpprefix}*

