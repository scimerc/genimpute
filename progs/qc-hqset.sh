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
  * While sex imputation is improved
    * Compile list of sex hq-variants
    * Compile list of non-sex hq-variants from input file
    * Extract all hq-variants from input file and make hq plink set
    * LD-prune hq-variants genotype data
    * Impute sex
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
[ -s "${genblacklist}" ] && declare -r regexcludeflag="--exclude range ${genblacklist}" || {
  printf "> error: file '%s' empty or not found.\n" "${genblacklist}" >&2;
  exit 1;
}

printf "> extracting high quality set..\n"

cp "${opt_inprefix}.bed" "${tmpprefix}_draft.bed"
cp "${opt_inprefix}.bim" "${tmpprefix}_draft.bim"
cp "${opt_inprefix}.fam" "${tmpprefix}_draft.fam"

# initialize hq variant set
> "${opt_outprefix}.mrk"

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
declare k=1
declare tmpindex=0
declare -a tmp_varmiss=( 0.01 0.05 0.1 )
declare Ntot
declare Nmin
declare Nhqi
Ntot=$( wc -l "${tmpprefix}_draft.bim" | tabulate | cut -f 1 )
Nmin=$(( Ntot / 3 ))
Nhqi=0
readonly Ntot
readonly Nmin
declare freqflag=''
declare hwetmpflag=''
declare imputesex=true
declare indcount=$( cat "${tmpprefix}_draft.fam" | wc -l )
if [ ${indcount} -lt ${cfg_minindcount} ] ; then
  [ -s "${opt_refprefix}.frq" ] && freqflag="--read-freq ${opt_refprefix}.frq" || {
    printf "> warning: too few individuals for accurate allele frequency estimates.\n"
    printf "> skipping sex imputation..\n"
    imputesex=false
  }
fi
declare qcdiff='qcdiff' # just some non-empty initialization string

while [ $Nhqi -lt $Nmin -a $tmpindex -le 2 -a "${qcdiff}" != "" ] ; do
  printf "> round $k..\n"
  # hq sex chromosomes
  xindcount=$( awk '$5 == 1 || $5 == 2' "${tmpprefix}_draft.fam" | wc -l )
  [ ${xindcount} -le ${cfg_minindcount} ] && hwetmpflag=''
  echo "  ${plinkexec} --allow-extra-chr --chr 23,24
               --bfile ${tmpprefix}_draft ${keepflag} ${extractflag} ${regexcludeflag}
               --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag}
               --make-just-bim
               --out ${tmpprefix}_sex" | printlog 2
  ${plinkexec} --allow-extra-chr --chr 23,24 \
               --bfile ${tmpprefix}_draft ${keepflag} ${extractflag} ${regexcludeflag} \
               --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag} \
               --make-just-bim \
               --out "${tmpprefix}_sex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # hq non-sex chromosomes
  hwetmpflag="--hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag}"
  echo "  ${plinkexec} --allow-extra-chr --not-chr 23,24
               --bfile ${tmpprefix}_draft ${keepflag} ${extractflag} ${regexcludeflag}
               --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag}
               --make-just-bim
               --out ${tmpprefix}_nonsex" | printlog 2
  ${plinkexec} --allow-extra-chr --not-chr 23,24 \
               --bfile ${tmpprefix}_draft ${keepflag} ${extractflag} ${regexcludeflag} \
               --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag} \
               --make-just-bim \
               --out "${tmpprefix}_nonsex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  cut -f 2 "${tmpprefix}_"*sex.bim | sort -u > "${tmpprefix}.mrk"
  # LD-prune hq variants
  echo "  ${plinkexec} --allow-extra-chr
               --bfile ${tmpprefix}_draft ${cfg_pruneflags}
               --extract "${tmpprefix}.mrk"
               --out ${tmpprefix}_ldp" | printlog 2
  ${plinkexec} --allow-extra-chr \
               --bfile "${tmpprefix}_draft" ${cfg_pruneflags} \
               --extract "${tmpprefix}.mrk" \
               --out "${tmpprefix}_ldp" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  sort -u "${tmpprefix}_ldp.prune.in" > "${tmpprefix}.imrk"
  sort -u -t $'\t' -k 2,2 "${tmpprefix}_"*sex.bim \
    | join -t $'\t' -2 2 "${tmpprefix}.imrk" - \
    > "${tmpprefix}_ldp.bim"
  tmpvarfile=''
  Nhq=$( sort -u "${tmpprefix}"*.mrk | wc -l )
  Nhqi=$( sort -u "${tmpprefix}"*.imrk | wc -l )
  NhqX=$( get_xvar_count "${tmpprefix}_"*sex.bim )
  NhqiX=$( get_xvar_count "${tmpprefix}_ldp.bim" )
  [ ${NhqX} -gt ${cfg_minvarcount} ] && tmpvarfile="${tmpprefix}.mrk"
  [ ${NhqiX} -gt ${cfg_minvarcount} ] && tmpvarfile="${tmpprefix}.imrk"
  # initialize sexcheck file
  awk 'BEGIN{
    OFS="\t"
    print( "FID", "IID", "SEXCHECK" )
  } { print( $1, $2, "__NA__" ) }' \
    "${tmpprefix}_draft.fam" > "${tmpprefix}_isex.sexcheck"
  # if there are enough X chromosome variants impute sex
  ${imputesex} && [ "${tmpvarfile}" != "" ] || break
  echo "  ${plinkexec} --allow-extra-chr
               --bfile ${tmpprefix}_draft ${freqflag}
               --extract "${tmpvarfile}"
               --impute-sex ycount
               --set-hh-missing
               --make-bed
               --out ${tmpprefix}_isex" | printlog 2
  ${plinkexec} --allow-extra-chr \
               --bfile "${tmpprefix}_draft" ${freqflag} \
               --extract "${tmpvarfile}" \
               --impute-sex ycount \
               --set-hh-missing \
               --make-bed \
               --out "${tmpprefix}_isex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # update sex in the draft file
  ${plinkexec} --allow-extra-chr \
               --bfile "${tmpprefix}_draft" \
               --update-sex "${tmpprefix}_isex.fam" 3 \
               --make-bed \
               --out "${tmpprefix}_out" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
  qcdiff=$( {
    cmp <( sort -u "${tmpprefix}_out.fam" ) <( sort -u "${tmpprefix}_draft.fam" ) 2>&1
    cmp <( sort -u "${tmpprefix}.mrk" ) <( sort -u "${opt_outprefix}.mrk" ) 2>&1
    } || true 
  )
  rename "${tmpprefix}_out" "${tmpprefix}_draft" "${tmpprefix}_out."*
  sort -u "${tmpprefix}.mrk" > "${opt_outprefix}.mrk"
  [ "${qcdiff}" == "" ] && tmpindex=$((tmpindex+1))
  k=$((k+1))
done
[ $Nhqi -gt $cfg_minvarcount ] || {
  printf "> error: not enough variants (%d) left in high quality set." $Nhqi >&2;
  exit 1;
}
# extract hq variants
${plinkexec} --allow-extra-chr \
             --bfile "${tmpprefix}_draft" \
             --extract "${tmpprefix}.mrk" \
             --set-hh-missing \
             --make-bed \
             --out "${tmpprefix}_out" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# extract LD-pruned hq variants
${plinkexec} --allow-extra-chr \
             --bfile "${tmpprefix}_draft" \
             --extract "${tmpprefix}.imrk" \
             --set-hh-missing \
             --make-bed \
             --out "${tmpprefix}_outi" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out"*.fam
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out"*.bim

mv "${tmpprefix}_isex.sexcheck" "${tmpprefix}_out.sexcheck"
rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out"*

rm -f "${tmpprefix}"*

