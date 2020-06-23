#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare -r tmpprefix="${opt_outprefix}_tmp"

declare -r cfg_fmax=$( cfgvar_get fmax )
declare -r cfg_fmin=$( cfgvar_get fmin )
declare -r cfg_freq=$( cfgvar_get freqhq )
declare -r cfg_genomeblacklist=$( cfgvar_get genomeblacklist )
declare -r cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_hweflag=$( cfgvar_get hweflag )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_pruneflags=$( cfgvar_get pruneflags )
declare -r cfg_samplemiss=$( cfgvar_get samplemiss )
declare -r cfg_uid=$( cfgvar_get uid )

#-------------------------------------------------------------------------------

# input: merged plink set
# output: hq plink set (with imputed sex, if possible)

echo -e "==== Variant HQC ====\n" | printlog 0

printf "\
  * While sex imputation improves
    * Compile list of sex hq-variants
    * Compile list of non-sex hq-variants from input file
    * Extract all hq-variants from input file and make hq plink set
    * LD-prune hq-variants genotype data
    * Impute sex (eventually)
\n" | printlog 0

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "> '%s' already exists. skipping hq..\n\n" "${opt_outprefix}.bed"
  exit 0
fi

if ls "${tmpprefix}"* > /dev/null 2>&1; then
  printf "> temporary files '%s*' found. please remove them before re-run.\n\n" "${tmpprefix}" >&2
  exit 1 
fi

#-------------------------------------------------------------------------------

declare keepflag=''
# set keep flag if a list of unrelated individuals exists
if [ -f "${opt_outprefixbase}/.i/qc/e_indqc.ids" ] ; then
  keepflag="--keep ${opt_outprefixbase}/.i/qc/e_indqc.ids"
fi
declare extractflag=''
# set extract flag if a previous list of variants exists
if [ -f "${opt_outprefixbase}/.i/qc/g_finqc.mrk" ] ; then
  extractflag="--extract ${opt_outprefixbase}/.i/qc/g_finqc.mrk"
fi

if [[ "${cfg_genomeblacklist}" == /* ]] ; then
  declare -r genblacklist="${cfg_genomeblacklist}"
else
  declare -r genblacklist="${BASEDIR}/lib/data/${cfg_genomeblacklist}"
fi
# check if exclude file exists and is not empty
[ -s "${genblacklist}" ] && declare -r regexcludeflag="--exclude range ${genblacklist}" || {
  printf "> error: file '%s' empty or not found.\n\n" "${genblacklist}" >&2;
  exit 1;
}

printf "> extracting high quality set..\n\n"

cp "${opt_inprefix}.bed" "${tmpprefix}_draft.bed"
cp "${opt_inprefix}.bim" "${tmpprefix}_draft.bim"
cp "${opt_inprefix}.fam" "${tmpprefix}_draft.fam"

# initialize hq variant sets
> "${tmpprefix}_draft.mrk"
> "${tmpprefix}_draft.imrk"

# if requested, merge with reference
if [ ! -z "${1:+x}" ] ; then
  ${plinkexec} --allow-extra-chr --bfile "${opt_refprefix}" \
               --extract <( cut -f 2 "${opt_inprefix}.bim" ) \
               --make-bed \
               --out "${tmpprefix}_ref" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  ${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" \
               --bmerge "${tmpprefix}_ref" \
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
    printf "> suppressing sex imputation..\n\n"
    imputesex=false
  }
fi
declare qcdiff='qcdiff' # just some non-empty initialization string to start the while clause

while [ $Nhqi -lt $Nmin -a $tmpindex -le ${#tmp_varmiss} -a "${qcdiff}" != "" ] ; do
  qcdiff=''
  printf "> round $k..\n" | printlog 1
  # hq sex chromosomes
  Nsexind=$( awk '$5 == 1 || $5 == 2' "${tmpprefix}_draft.fam" | wc -l )
  Nsexvar=$( awk '$1 == 23 || $1 == 24' "${tmpprefix}_draft.bim" | wc -l )
  # deactivate HWE tests for sex chromosomes if too few individuals have pre-determined sex
  [ ${Nsexind} -le ${cfg_minindcount} ] && hwetmpflag=''
  # extract set of high quality sex chromosome variants if enough are there to begin with
  if [ ${Nsexvar} -gt 0 ] ; then
    echo -e "  ${plinkexec##*/} --allow-extra-chr --chr 23,24
                 --bfile ${tmpprefix}_draft ${keepflag} ${extractflag} ${regexcludeflag}
                 --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag}
                 --mind ${cfg_samplemiss}
                 --make-just-bim
                 --out ${tmpprefix}_sex\n" | printlog 2
    ${plinkexec} --allow-extra-chr --chr 23,24 \
                 --bfile "${tmpprefix}_draft" ${keepflag} ${extractflag} ${regexcludeflag} \
                 --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag} \
                 --mind ${cfg_samplemiss} \
                 --make-just-bim \
                 --out "${tmpprefix}_sex" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
  fi
  # hq non-sex chromosomes
  hwetmpflag="--hwe 1.E-${cfg_hweneglogp_ctrl} ${cfg_hweflag}"
  # extract set of high quality non-sex chromosome variants
  echo -e "  ${plinkexec##*/} --allow-extra-chr --not-chr 23,24
               --bfile ${tmpprefix}_draft ${keepflag} ${extractflag} ${regexcludeflag}
               --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag}
               --mind ${cfg_samplemiss}
               --make-just-bim
               --out ${tmpprefix}_nonsex\n" | printlog 2
  ${plinkexec} --allow-extra-chr --not-chr 23,24 \
               --bfile "${tmpprefix}_draft" ${keepflag} ${extractflag} ${regexcludeflag} \
               --geno ${tmp_varmiss[${tmpindex}]} ${hwetmpflag} --maf ${cfg_freq} ${freqflag} \
               --mind ${cfg_samplemiss} \
               --make-just-bim \
               --out "${tmpprefix}_nonsex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  cut -f 2 "${tmpprefix}_"*sex.bim | sort -u > "${tmpprefix}_devel.mrk"
  # LD-prune hq variants
  echo -e "  ${plinkexec##*/} --allow-extra-chr
               --bfile ${tmpprefix}_draft ${cfg_pruneflags}
               --extract ${tmpprefix}_devel.mrk
               --out ${tmpprefix}_ldp\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
               --bfile "${tmpprefix}_draft" ${cfg_pruneflags} \
               --extract "${tmpprefix}_devel.mrk" \
               --out "${tmpprefix}_ldp" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # ensure a prune.in file exists
  touch "${tmpprefix}_ldp.prune.in"
  sort -u "${tmpprefix}_ldp.prune.in" > "${tmpprefix}_devel.imrk"
  sort -u -t $'\t' -k 2,2 "${tmpprefix}_"*sex.bim \
    | join -t $'\t' -2 2 "${tmpprefix}_devel.imrk" - \
    > "${tmpprefix}_ldp.bim"
  tmpfflag=''
  tmpvarfile=''
  Nhq=$( sort -u "${tmpprefix}"*.mrk | wc -l )
  Nhqi=$( sort -u "${tmpprefix}"*.imrk | wc -l )
  NhqX=$( get_xvar_count "${tmpprefix}_"*sex.bim )
  NhqiX=$( get_xvar_count "${tmpprefix}_ldp.bim" )
  if [ ${NhqX} -gt ${cfg_minvarcount} ] ; then
    tmpvarfile="${tmpprefix}_devel.mrk"
  fi
  if [ ${NhqiX} -gt ${cfg_minvarcount} ] ; then
    tmpfflag="${cfg_fmax} ${cfg_fmin}"
    tmpvarfile="${tmpprefix}_devel.imrk"
  fi
  # initialize sexcheck file
  awk 'BEGIN{
    OFS="\t"
    print( "FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F_SEX" )
  } { print( $1, $2, "__NA__", "__NA__", "__NA__", "__NA__" ) }' \
    "${tmpprefix}_draft.fam" > "${tmpprefix}_isex.sexcheck"
  # if there are enough X chromosome variants impute sex and update, else copy
  if [ ${imputesex} -a "${tmpvarfile}" != "" ] ; then
    echo -e "  ${plinkexec##*/} --allow-extra-chr
                 --bfile ${tmpprefix}_draft ${freqflag}
                 --extract ${tmpvarfile}
                 --impute-sex ycount ${tmpfflag}
                 --set-hh-missing
                 --make-bed
                 --out ${tmpprefix}_isex\n" | printlog 2
    ${plinkexec} --allow-extra-chr \
                 --bfile "${tmpprefix}_draft" ${freqflag} \
                 --extract "${tmpvarfile}" \
                 --impute-sex ycount ${tmpfflag} \
                 --set-hh-missing \
                 --make-bed \
                 --out "${tmpprefix}_isex" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    # update sex in the draft file
    echo -e "  ${plinkexec##*/} --allow-extra-chr
                 --bfile ${tmpprefix}_draft
                 --update-sex ${tmpprefix}_isex.fam 3
                 --make-bed
                 --out ${tmpprefix}_devel\n" | printlog 2
    ${plinkexec} --allow-extra-chr \
                 --bfile "${tmpprefix}_draft" \
                 --update-sex "${tmpprefix}_isex.fam" 3 \
                 --make-bed \
                 --out "${tmpprefix}_devel" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_devel.fam"
    sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_devel.bim"
    qcdiff=$( {
      cmp <( sort -u "${tmpprefix}_devel.fam" ) <( sort -u "${tmpprefix}_draft.fam" ) 2>&1
      cmp <( sort -u "${tmpprefix}_devel.mrk" ) <( sort -u "${tmpprefix}_draft.mrk" ) 2>&1
      } || true )
  else
    cp "${tmpprefix}_draft.bed" "${tmpprefix}_devel.bed"
    cp "${tmpprefix}_draft.bim" "${tmpprefix}_devel.bim"
    cp "${tmpprefix}_draft.fam" "${tmpprefix}_devel.fam"
  fi
  # reset draft set to newly created devel set..
  cp "${tmpprefix}_devel.bed" "${tmpprefix}_draft.bed"
  cp "${tmpprefix}_devel.bim" "${tmpprefix}_draft.bim"
  cp "${tmpprefix}_devel.fam" "${tmpprefix}_draft.fam"
  # ..and copy the corresponding variant sets
  cp "${tmpprefix}_devel.mrk" "${tmpprefix}_draft.mrk"
  cp "${tmpprefix}_devel.imrk" "${tmpprefix}_draft.imrk"
  # if there aren't any changes increase tmpindex for an eventual additional round, otherwise
  # re-run with the same parameters, re-impute sex and re-check back here for new changes
  [ "${qcdiff}" == "" ] && tmpindex=$((tmpindex+1))
  # increase round counter
  k=$((k+1))
done
[ $Nhqi -gt $cfg_minvarcount ] || {
  printf "> not enough variants (%d) left in high quality set: aborting..\n\n" $Nhqi >&2;
  exit 1;
}
# extract hq variants
echo -e "  ${plinkexec##*/} --allow-extra-chr
             --bfile ${tmpprefix}_draft
             --extract ${tmpprefix}_draft.mrk
             --set-hh-missing
             --make-bed
             --out ${tmpprefix}_out\n" | printlog 2
${plinkexec} --allow-extra-chr \
             --bfile "${tmpprefix}_draft" \
             --extract "${tmpprefix}_draft.mrk" \
             --set-hh-missing \
             --make-bed \
             --out "${tmpprefix}_out" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# extract LD-pruned hq variants
echo -e "  ${plinkexec##*/} --allow-extra-chr
             --bfile ${tmpprefix}_draft
             --extract ${tmpprefix}_draft.imrk
             --set-hh-missing
             --make-bed
             --out ${tmpprefix}_outi\n" | printlog 2
${plinkexec} --allow-extra-chr \
             --bfile "${tmpprefix}_draft" \
             --extract "${tmpprefix}_draft.imrk" \
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

