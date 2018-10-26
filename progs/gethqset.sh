#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_regionblacklist=${BASEDIR}/data/highLD_b37.bed
declare -r cfg_varmiss=0.05
declare -r cfg_freqhq=0.2
declare -r cfg_hweneglogp_ctrl=4
declare -r cfg_hweflag='midp include-nonctrl'
declare -r cfg_minindcount=100
declare -r cfg_minvarcount=100
declare -r cfg_uid=000UID

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
# 3) extract all hq variants from input file and make hq plink set
# 4) LD-prune hq variants from 3)
# 5) impute sex once with all standard hq variants from 4
# 6) if sex could be imputed for enough individuals, then
#      impute it once again after HWE tests
# 7) update biography file with sex information


# check if exclude file exists and is not empty
[ -s "${cfg_regionblacklist}" ] || {
  printf "error: file '%s' empty or not found.\n" ${cfg_regionblacklist} >&2;
  exit 1;
}
declare -r excludeopt="--exclude range ${cfg_regionblacklist}"

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

tabulate() {
  sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;'
}

paste_sample_ids() {
  local -r infile="$1"
  tabulate < "${infile}" \
    | awk -F $'\t' -v uid=${cfg_uid} '{
        OFS="\t"
        if ( NR>1 ) uid=$1"_"$2
        printf( "%s", uid )
        for ( k=3; k<=NF; k++ )
          printf( "\t%s", $k )
        printf( "\n" )
      }' \
    | sort -t $'\t' -u -k 1,1
}

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

  declare -r xindcount=$( awk '$5 == 1 || $5 == 2' ${tmpprefix}_hq_LDpruned_isex.fam | wc -l )
  # if sex could be imputed for enough individuals impute it once again after HWE tests
  if (( xindcount > minindcount )) ; then
    plink --bfile ${tmpprefix}_hq_LDpruned_isex \
          --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hwflag} \
          --make-just-bim \
          --out ${tmpprefix}_hq_LDpruned_sexhwe \
          >> ${debuglogfn}

    # if there are enough X chromosome variants after HWE re-impute sex based on them
    if [ $( get_xvar_count ${tmpprefix}_hq_LDpruned_sexhwe.bim ) -gt $cfg_minvarcount ] ; then
      plink --bfile ${tmpprefix}_hq_LDpruned_isex \
            --extract <( cut -f2 ${tmpprefix}_hq_LDpruned_sexhwe.bim ) \
            --impute-sex \
            --make-bed \
            --out ${tmpprefix}_hq_LDpruned_isex_new \
            >> ${debuglogfn}
      # replace the original sex imputation files
      mv ${tmpprefix}_hq_LDpruned_isex_new.bed ${tmpprefix}_hq_LDpruned_isex.bed
      mv ${tmpprefix}_hq_LDpruned_isex_new.bim ${tmpprefix}_hq_LDpruned_isex.bim
      mv ${tmpprefix}_hq_LDpruned_isex_new.fam ${tmpprefix}_hq_LDpruned_isex.fam
      mv ${tmpprefix}_hq_LDpruned_isex_new.sexcheck ${tmpprefix}_hq_LDpruned_isex.sexcheck
    fi
  fi
  {
    paste_sample_ids ${tmpprefix}_hq_LDpruned_isex.sexcheck \
      | join -t $'\t' ${opt_biofile} - \
      | tee ${tmpprefix}.0.bio
    TNF=$( wc -l ${tmpprefix}.0.bio | tabulate | cut -f 1 )
    paste_sample_ids ${tmpprefix}_hq_LDpruned_isex.sexcheck \
      | join -t $'\t' -v1 ${opt_biofile} - \
      | awk -F $'\t' -v TNF=${TNF} '{
        OFS="\t"
        printf($0)
        for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
        printf("\n")
      }'
  } | sort -u -k 1,1 > ${tmpprefix}.1.bio
  cp ${tmpprefix}.1.bio ${opt_biofile}
fi

mv ${tmpprefix}_hq_LDpruned_isex.bed ${opt_outprefix}.bed
mv ${tmpprefix}_hq_LDpruned_isex.bim ${opt_outprefix}.bim
mv ${tmpprefix}_hq_LDpruned_isex.fam ${opt_outprefix}.fam

# rm ${tmpprefix}*

