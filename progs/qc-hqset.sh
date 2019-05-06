#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare -r tmpprefix="${opt_outprefix}_tmp"

declare -r cfg_freq=$( cfgvar_get freqhq )
declare -r cfg_genomeblacklist=$( cfgvar_get genomeblacklist )
declare -r cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_pruneflags=$( cfgvar_get pruneflags )
declare -r cfg_uid=$( cfgvar_get uid )

#-------------------------------------------------------------------------------

# input: merged plink set
# output: hq plink set (with imputed sex, if possible)

echo -e "==== Variant HQC ====\n" | printlog 1

printf "\
  * Compile list of sex hq-variants
  * Compile list of non-sex hq-variants from input file
  * Extract all hq-variants from input file and make hq plink set
  * LD-prune hq-variants genotype data
  * Impute sex once with all standard hq-variants
  * If sex could be imputed for enough individuals, impute it once again after \
    sex-chromosome HWE tests
\n" | printlog 1

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "> '%s' already exists. skipping hq..\n" "${opt_outprefix}.bed"
  exit 0
fi

if ls "${tmpprefix}"* > /dev/null 2>&1; then
  printf "> temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1 
fi

#-------------------------------------------------------------------------------

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

declare -r genblacklist="${BASEDIR}/lib/data/${cfg_genomeblacklist}"
# check if exclude file exists and is not empty
[ -s "${genblacklist}" ] || {
  printf "> error: file '%s' empty or not found.\n" "${genblacklist}" >&2;
  exit 1;
}
declare -r regexcludeflag="--exclude range ${genblacklist}"

printf "> extracting high quality set..\n"

cp "${opt_inprefix}.bed" "${tmpprefix}_draft.bed"
cp "${opt_inprefix}.bim" "${tmpprefix}_draft.bim"
cp "${opt_inprefix}.fam" "${tmpprefix}_draft.fam"

if [ ! -z "${1:+x}" ] ; then
  # if requested merge with reference
  ${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" \
               --bmerge "${opt_refprefix}" \
               --out "${tmpprefix}_draft" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
fi
declare tmpindex=0
declare -a tmp_varmiss=( 0.01 0.05 0.1 )
declare Ntot
declare Nmin
declare Nhq
Ntot=$( wc -l "${tmpprefix}_draft.bim" | tabulate | cut -f 1 )
Nmin=$(( Ntot / 3 ))
Nhq=0
# enable progressively milder coverage requirements to accomodate smaller batches
while [ $Nhq -lt $Nmin -a $tmpindex -le 2 ] ; do
  # get non sex hq-variants from input file
  echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_draft ${keepflag}
               --not-chr 23,24 ${regexcludeflag} ${extractflag}
               --geno ${tmp_varmiss[${tmpindex}]}
               --maf ${cfg_freq}
               --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag}
               --make-just-bim
               --out ${tmpprefix}_nonsex" | printlog 2
  ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_draft" ${keepflag} \
               --not-chr 23,24 ${regexcludeflag} ${extractflag} \
               --geno ${tmp_varmiss[${tmpindex}]} \
               --maf ${cfg_freq} \
               --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag} \
               --make-just-bim \
               --out "${tmpprefix}_nonsex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  if [ $( get_xvar_count "${tmpprefix}_draft.bim" ) -ge ${cfg_minvarcount} ] ; then
    # get sex hq-variants from input file
    echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_draft ${keepflag}
                 --chr 23,24 ${regexcludeflag} ${extractflag}
                 --geno ${tmp_varmiss[${tmpindex}]}
                 --maf ${cfg_freq}
                 --make-just-bim
                 --out ${tmpprefix}_sex" | printlog 2
    ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_draft" ${keepflag} \
                 --chr 23,24 ${regexcludeflag} ${extractflag} \
                 --geno ${tmp_varmiss[${tmpindex}]} \
                 --maf ${cfg_freq} \
                 --make-just-bim \
                 --out "${tmpprefix}_sex" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
  fi
  Nhq=$( sort -u -k 2,2 "${tmpprefix}_"*sex.bim | wc -l )
  tmpindex=$(( tmpindex + 1 ))
done
# check if we have anything of high quality
[ $Nhq -ge $cfg_minvarcount ] || {
  printf "> error: not enough variants left in high quality set." >&2;
  exit 1;
}
# extract all hq variants from input file and make hq plink set
cut -f 2 "${tmpprefix}_"*sex.bim | sort -u > "${tmpprefix}_hq.mrk"
${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_draft" \
             --extract "${tmpprefix}_hq.mrk" \
             --make-bed \
             --out "${tmpprefix}_hq" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# LD-prune hq variants
echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_hq ${keepflag} ${cfg_pruneflags}
             --out ${tmpprefix}_hq_LD" | printlog 2
${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_hq" ${keepflag} ${cfg_pruneflags} \
             --out "${tmpprefix}_hq_LD" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
# extract LD-pruned hq variants from hq plink set
${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_hq" \
             --extract "${tmpprefix}_hq_LD.prune.in" \
             --make-bed \
             --out "${tmpprefix}_hq_LDpruned" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# if there are enough X chromosome variants impute sex based on them
if [ $( get_xvar_count "${tmpprefix}_hq_LDpruned.bim" ) -ge ${cfg_minvarcount} ] ; then
  declare freqflag=''
  declare -r indcount=$( cat "${tmpprefix}_hq_LDpruned.fam" | wc -l )
  if [ ${indcount} -lt ${cfg_minindcount} ] ; then
    printf "> warning: insufficient individuals for accurate allele frequency estimates.\n"
    [ -s "${opt_refprefix}.frq" ] && freqflag="--read-freq ${opt_refprefix}.frq"
  fi
  # impute sex once with all standard high quality variants
  echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_hq_LDpruned ${freqflag}
               --impute-sex ycount
               --make-bed
               --out ${tmpprefix}_hq_LDpruned_isex" | printlog 2
  ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_hq_LDpruned" ${freqflag} \
               --impute-sex ycount \
               --make-bed \
               --out "${tmpprefix}_hq_LDpruned_isex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  rename _hq_LDpruned_isex _out "${tmpprefix}_hq_LDpruned_isex".*
  declare -r xindcount=$( awk '$5 == 1 || $5 == 2' "${tmpprefix}_out.fam" | wc -l )
  # if sex could be imputed for enough individuals impute it once again after HWE tests
  if [ ${xindcount} -ge ${cfg_minindcount} ] ; then
    echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_out
                 --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag}
                 --make-just-bim
                 --out ${tmpprefix}_sexhwe" | printlog 2
    ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_out" \
                 --hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag} \
                 --make-just-bim \
                 --out "${tmpprefix}_sexhwe" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    sed -r 's/[ \t]+/\t/g' "${tmpprefix}_sexhwe.bim" \
      | cut -f 2 > "${tmpprefix}_sexhwe.mrk"
    # if there are enough X chromosome variants after HWE re-impute sex based on them
    if [ $( get_xvar_count "${tmpprefix}_sexhwe.bim" ) -ge ${cfg_minvarcount} ] ; then
      echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_hq_LDpruned ${freqflag}
                   --extract ${tmpprefix}_sexhwe.mrk
                   --impute-sex ycount
                   --make-bed
                   --out ${tmpprefix}_hq_LDpruned_isexhwe" | printlog 2
      ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_hq_LDpruned" ${freqflag} \
                   --extract <( cut -f 2 "${tmpprefix}_sexhwe.bim" ) \
                   --impute-sex ycount \
                   --make-bed \
                   --out "${tmpprefix}_hq_LDpruned_isexhwe" \
                   2> >( tee "${tmpprefix}.err" ) | printlog 3
      if [ $? -ne 0 ] ; then
        cat "${tmpprefix}.err"
      fi
      # replace the original sex imputation files
      rename _hq_LDpruned_isexhwe _out "${tmpprefix}_hq_LDpruned_isexhwe".*
    fi
  fi
else
  # if sex could not be imputed use LD-pruned set
  rename _hq_LDpruned _out "${tmpprefix}_hq_LDpruned".*
  # write a dummy sexcheck file for later use
  awk 'BEGIN{
    OFS="\t"
    print( "FID", "IID", "SEXCHECK" )
  } { print( $1, $2, "__NA__" ) }' \
  "${tmpprefix}_out.fam" > "${tmpprefix}_out.sexcheck"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"

rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out".*

rm -f "${tmpprefix}"*

