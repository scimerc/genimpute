#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r  debuglogfn=${opt_outprefix}_tmp_debug.log

if [ -f "${opt_hqprefix}.genome.gz" -a -f "${opt_hqprefix}.eigenvec" ] ; then
  printf "skipping PCA..\n"
  exit 0
fi

# input: high quality plink set
# output: gzipped GRM and a matrix of eigenvector weights


plink --bfile ${opt_hqprefix} --genome gz --out ${opt_hqprefix} >> ${debuglogfn}

plink --bfile ${opt_hqprefix} \
      --cluster \
      --read-genome ${opt_hqprefix}.genome.gz \
      --pca header tabs var-wts \
      --out ${opt_hqprefix} \
      >> ${debuglogfn}

