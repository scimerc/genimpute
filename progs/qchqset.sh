#!/usr/bin/env bash

# exit on error
trap 'printf "=== error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_varmiss=$( cfgvar_get varmiss )
declare -r cfg_freqhq=$( cfgvar_get freqhq )
declare -r cfg_genomeblacklist=$( cfgvar_get genomeblacklist )
declare -r cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_uid=$( cfgvar_get uid )

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping hq step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
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


declare -r regionblacklist=${BASEDIR}/lib/data/${cfg_genomeblacklist}
# check if exclude file exists and is not empty
[ -s "${regionblacklist}" ] || {
  printf "error: file '%s' empty or not found.\n" ${regionblacklist} >&2;
  exit 1;
}
declare -r excludeopt="--exclude range ${regionblacklist}"

# get sex hq-variants from input file
plink --bfile ${opt_inprefix} \
      --not-chr 23,24 ${excludeopt} \
      --geno ${cfg_varmiss} \
      --maf ${cfg_freqhq} \
      --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag} \
      --make-just-bim \
      --out ${tmpprefix}_nonsex \
      >> ${debuglogfn}

# get non-sex hq-variants from input file
plink --bfile ${opt_inprefix} \
      --chr 23,24 ${excludeopt} \
      --geno ${cfg_varmiss} \
      --maf ${cfg_freqhq} \
      --make-just-bim \
      --out ${tmpprefix}_sex \
      >> ${debuglogfn}


# check if we have anything of high quality
[ -s "${tmpprefix}_nonsex.bim" -o -s "${tmpprefix}_sex.bim" ] || {
  printf "error: no variants left in high quality set." >&2;
  exit 1;
}

# extract all hq variants from input file and make hq plink set
cut -f 2 ${tmpprefix}_*sex.bim | sort -u > ${tmpprefix}_hq.mrk
plink --bfile ${opt_inprefix} \
      --extract ${tmpprefix}_hq.mrk \
      --make-bed \
      --out ${tmpprefix}_hq \
      >> ${debuglogfn}

# LD-prune hq variants
plink --bfile ${tmpprefix}_hq \
      --indep-pairphase 500 5 0.2 \
      --out ${tmpprefix}_hq_LD \
      >> ${debuglogfn}

# extract LD-pruned hq variants from hq plink set
plink --bfile ${tmpprefix}_hq \
      --extract ${tmpprefix}_hq_LD.prune.in \
      --make-bed \
      --out ${tmpprefix}_hq_LDpruned \
      >> ${debuglogfn}

get_xvar_count() {
  awk '$1 == 23' $1 | wc -l
}

# if there are enough X chromosome variants impute sex based on them
if [ $( get_xvar_count ${tmpprefix}_hq_LDpruned.bim ) -gt $cfg_minvarcount ] ; then
  # impute sex once with all standard high quality variants
  plink --bfile ${tmpprefix}_hq_LDpruned \
        --impute-sex \
        --make-bed \
        --out ${tmpprefix}_hq_LDpruned_isex \
        >> ${debuglogfn}
  mv ${tmpprefix}_hq_LDpruned_isex.bed ${tmpprefix}_out.bed
  mv ${tmpprefix}_hq_LDpruned_isex.bim ${tmpprefix}_out.bim
  mv ${tmpprefix}_hq_LDpruned_isex.fam ${tmpprefix}_out.fam
  mv ${tmpprefix}_hq_LDpruned_isex.sexcheck ${tmpprefix}_out.sexcheck

  declare -r xindcount=$( awk '$5 == 1 || $5 == 2' ${tmpprefix}_out.fam | wc -l )
  # if sex could be imputed for enough individuals impute it once again after HWE tests
  if (( xindcount > minindcount )) ; then
    plink --bfile ${tmpprefix}_out \
          --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hwflag} \
          --make-just-bim \
          --out ${tmpprefix}_sexhwe \
          >> ${debuglogfn}

    # if there are enough X chromosome variants after HWE re-impute sex based on them
    if [ $( get_xvar_count ${tmpprefix}_sexhwe.bim ) -gt $cfg_minvarcount ] ; then
      plink --bfile ${tmpprefix}_out \
            --extract <( cut -f 2 ${tmpprefix}_sexhwe.bim ) \
            --impute-sex \
            --make-bed \
            --out ${tmpprefix}_hq_LDpruned_isex_new \
            >> ${debuglogfn}
      # replace the original sex imputation files
      mv ${tmpprefix}_hq_LDpruned_isex_new.bed ${tmpprefix}_out.bed
      mv ${tmpprefix}_hq_LDpruned_isex_new.bim ${tmpprefix}_out.bim
      mv ${tmpprefix}_hq_LDpruned_isex_new.fam ${tmpprefix}_out.fam
      mv ${tmpprefix}_hq_LDpruned_isex_new.sexcheck ${tmpprefix}_out.sexcheck
    fi
  fi
fi

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
mv ${tmpprefix}_out.sexcheck ${opt_outprefix}.sexcheck

rm ${tmpprefix}*

