#!/usr/bin/env bash

set -Eeou pipefail

declare -r  tmpprefix=${opt_outprefix}_tmp
declare -r  debuglogfn=${tmpprefix}_debug.log
declare -ra batchfiles=( ${opt_inputfiles} )

#-------------------------------------------------------------------------------

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping merge..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

#-------------------------------------------------------------------------------

# init global blacklist
# while true
  # for every batch
    # purge batchfile of global blacklist
  # attempt to merge resulting batchfiles
  # stop on success
  # if blacklist is not empty, tell francesco
  #   (if we run more than 2 times something might be weird)
  # append to global blacklist (derived from errors)


echo "merging batches.."
declare -r varblacklist=${tmpprefix}.exclude
declare -r mismatchlist=${tmpprefix}.mismatch
declare -r batchlist=${tmpprefix}.batchlist
while true ; do
  > ${batchlist}
  for i in ${!batchfiles[@]} ; do
    declare batchcode=$( get_unique_filename_from_path ${batchfiles[$i]} )
    declare b_inprefix=${opt_inprefix}_batch_${batchcode}
    declare b_outprefix=${opt_outprefix}_batch_${batchcode}
    declare plinkflag=''
    if [ -s "${varblacklist}" ] ; then
      plinkflag="--exclude ${varblacklist}" 
    fi
    plink --bfile ${b_inprefix} ${plinkflag} \
          --make-bed \
          --out ${b_outprefix} \
          2>&1 >> ${debuglogfn} \
          | tee -a ${debuglogfn}
    echo "${b_outprefix}" >> ${batchlist}
    unset batchcode
    unset inprefix
    unset outprefix
  done
  unset plinkflag
  # NOTE: if only one batch is present plink throws a warning
  plink --merge-list ${batchlist} \
        --out ${tmpprefix}_out \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  # extract plink's warnings about chromosome and position clashes from plink's log and add the
  # corresponding variant names to plink's own missnp file.
  egrep '^Warning: Multiple' ${tmpprefix}_out.log \
    | cut -d ' ' -f 7 \
    | tr -d "'." \
    | sort -u \
    >> ${mismatchlist}
  if [[ -s "${mismatchlist}" ]] ; then
    sort -u ${mismatchlist} >> ${tmpprefix}_out.missnp
  fi
  # are we done?
  if [ ! -f "${tmpprefix}.missnp" ] ; then
    break
  fi
  # no! prepare to repeat loop
  sort -u ${tmpprefix}.missnp > ${varblacklist}
  echo "repeating merging attempt.."
  rm ${tmpprefix}.missnp
done
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.fam ${opt_outprefix}.fam

rm ${tmpprefix}*

