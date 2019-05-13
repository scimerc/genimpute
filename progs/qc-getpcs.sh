#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare -r tmpprefix="${opt_outprefix}_tmp"

declare -r cfg_uid=$( cfgvar_get uid )

#-------------------------------------------------------------------------------

# input: high quality plink set
# output: gzipped GRM and a matrix of eigenvector weights

echo -e "==== PCA etc. ====\n" | printlog 1

printf "\
  * Compute genetic principal components
  * Compute final individual coverage statistics
\n" | printlog 1

if [ -f "${opt_outprefix}.eigenvec" -a -f "${opt_outprefix}.eigenvec.var" ] ; then
  printf "> '%s' found. skipping PCA..\n" "${opt_outprefix}.eigenvec"
  exit 0
fi

printf "> computing genetic principal components..\n"

${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}" \
             --genome gz \
             --out "${tmpprefix}_draft" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
mv "${tmpprefix}_draft.genome.gz" "${opt_outprefix}.genome.gz"
${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}" \
             --cluster \
             --read-genome "${opt_outprefix}.genome.gz" \
             --pca header tabs var-wts \
             --out "${tmpprefix}_draft" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
rename "${tmpprefix}_draft" "${opt_outprefix}" "${tmpprefix}_draft.eigen"*
${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" \
             --missing \
             --out "${tmpprefix}_draft" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
mv ${tmpprefix}_draft.imiss ${opt_outprefix}.imiss
mv ${tmpprefix}_draft.lmiss ${opt_outprefix}.lmiss

# update biography file with genetic PCs
{
  {
    printf "%s\t" ${cfg_uid}
    head -n 1 "${opt_outprefix}.eigenvec" | tabulate | cut -f 3-
  } | join -t $'\t'     "${opt_biofile}" - \
    | tee "${tmpprefix}.0.bio"
  # count number of fields in the merged file
  TNF=$( cat "${tmpprefix}.0.bio" | wc -w )
  # merge information for existing individuals
  attach_uids "${opt_outprefix}.eigenvec" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 | join -t $'\t'     "${opt_biofile}" - \
  # add non-existing individuals and pad the extra fields with NAs
  attach_uids "${opt_outprefix}.eigenvec" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.1.bio"
mv "${tmpprefix}.1.bio" "${opt_biofile}"

# update biography file with missingness statistics
{
  {
    printf "%s\t" ${cfg_uid}
    head -n 1 "${opt_outprefix}.imiss" | tabulate | cut -f 3-
  } | join -t $'\t'     "${opt_biofile}" - \
    | tee "${tmpprefix}.0.bio"
  # count number of fields in the merged file
  TNF=$( cat "${tmpprefix}.0.bio" | wc -w )
  # merge information for existing individuals
  attach_uids "${opt_outprefix}.imiss" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 | join -t $'\t'     "${opt_biofile}" - \
  # add non-existing individuals and pad the extra fields with NAs
  attach_uids "${opt_outprefix}.imiss" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.3.bio"
mv "${tmpprefix}.3.bio" "${opt_biofile}"

rm "${tmpprefix}"*

