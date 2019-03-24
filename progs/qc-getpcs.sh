#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_uid=$( cfgvar_get uid )

#-------------------------------------------------------------------------------

# input: high quality plink set
# output: gzipped GRM and a matrix of eigenvector weights

printf "\
  * Compute genetic principal components
  * Compute final individual coverage statistics
" | printlog 1

if [ -f "${opt_hqprefix}.eigenvec" -a -f "${opt_hqprefix}.eigenvec.var" ] ; then
  printf "'%s' found. skipping PCA..\n" "${opt_hqprefix}.eigenvec"
  exit 0
fi

${plinkexec} --bfile ${opt_hqprefix} \
             --genome gz \
             --out ${opt_hqprefix} \
             2>&1 | printlog 2
${plinkexec} --bfile ${opt_hqprefix} \
             --cluster \
             --read-genome ${opt_hqprefix}.genome.gz \
             --pca header tabs var-wts \
             --out ${opt_hqprefix} \
             2>&1 | printlog 2
${plinkexec} --bfile ${opt_inprefix} \
             --missing \
             --out ${opt_outprefix} \
             2>&1 | printlog 2

# update biography file with genetic PCs
{
  attach_uids ${opt_hqprefix}.eigenvec \
    | join -t $'\t'     ${opt_biofile} - \
    | tee ${tmpprefix}.0.bio
  TNF=$( head -n 1 ${tmpprefix}.0.bio | wc -w )
  attach_uids ${opt_hqprefix}.eigenvec \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.1.bio
mv ${tmpprefix}.1.bio ${opt_biofile}

# update biography file with missingness statistics
{
  attach_uids ${opt_outprefix}.imiss \
    | join -t $'\t' ${opt_biofile} - \
    | tee ${tmpprefix}.2.bio
  TNF=$( head -n 1 ${tmpprefix}.2.bio | wc -w )
  attach_uids ${opt_outprefix}.imiss \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.3.bio
mv ${tmpprefix}.3.bio ${opt_biofile}

rm ${tmpprefix}*

