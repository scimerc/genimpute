#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r  tmpprefix="${opt_outprefix}_tmp"
declare -r  debuglogfn="${tmpprefix}_debug.log"
declare -ra batchfiles=( ${opt_inputfiles} )

#-------------------------------------------------------------------------------

# input: aligned plink file sets
# output: merged plink file set, eventually purged of conflicts

printf "\
  * Initialize global blacklist
  * While true
    * For every batch
      * Purge batchfile of global blacklist (if any)
    * Attempt merging resulting batchfiles
    * Stop on success
    * Abort procedure if blacklist is not empty \
      (if we run more than 2 times something might be weird)
    * Append mismatching variants to global blacklist (derived from errors)
" | printlog 1

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "'%s' already exists. skipping merge..\n" "${opt_outprefix}.bed"
  exit 0
fi

if ls "${tmpprefix}"* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

#-------------------------------------------------------------------------------

printf "merging batches..\n"

declare -r varblacklist="${tmpprefix}.exclude"
declare -r mismatchlist="${tmpprefix}.mismatch"
declare -r batchlist="${tmpprefix}.batchlist"
while true ; do
  > "${batchlist}"
  for i in ${!batchfiles[@]} ; do
    declare batchcode=$( get_unique_filename_from_path ${batchfiles[$i]} )
    declare b_inprefix="${opt_inprefix}_batch_${batchcode}"
    declare b_outprefix="${opt_outprefix}_batch_${batchcode}"
    declare plinkflag=''
    if [ -s "${varblacklist}" ] ; then
      plinkflag="--exclude ${varblacklist}" 
    fi
    ${plinkexec} \
          --bfile "${b_inprefix}" ${plinkflag} \
          --make-bed \
          --out "${b_outprefix}" \
          2> >( tee "${tmpprefix}.err" ) | printlog 2
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    printf "${b_outprefix}\n" >> "${batchlist}"
    unset batchcode
    unset inprefix
    unset outprefix
  done
  unset plinkflag
  # NOTE: if only one batch is present plink throws a warning
  ${plinkexec} \
        --merge-list "${batchlist}" \
        --out "${tmpprefix}_out" \
        2> >( tee "${tmpprefix}.err" ) | printlog 2
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # extract plink's warnings about chromosome and position clashes from plink's log and add the
  # corresponding variant names to plink's own missnp file.
  grep '^Warning: Multiple' "${tmpprefix}_out.log" \
    | cut -d ' ' -f 7 \
    | tr -d "'." \
    | sort -u \
    >> "${mismatchlist}" || true
  if [[ -s "${mismatchlist}" ]] ; then
    sort -u "${mismatchlist}" >> "${tmpprefix}_out.missnp"
  fi
  # are we done?
  if [ ! -f "${tmpprefix}.missnp" ] ; then
    break
  fi
  # no! prepare to repeat loop
  sort -u "${tmpprefix}.missnp" > "${varblacklist}"
  printf "repeating merging attempt..\n"
  rm "${tmpprefix}.missnp"
done
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"

rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out".*

rm "${tmpprefix}"*

