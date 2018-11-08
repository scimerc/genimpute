#!/usr/bin/env bash

# exit on error
trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_uid=$( cfgvar_get uid )

if [ -f "${opt_hqprefix}.eigenvec" -a -f "${opt_hqprefix}.eigenvec.var" ] ; then
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

# update biography file with genetic PCs
{
  paste_sample_ids ${opt_hqprefix}.eigenvec \
    | join -t $'\t' ${opt_biofile} - \
    | tee ${tmpprefix}.0.bio
  TNF=$( wc -l ${tmpprefix}.0.bio | tabulate | cut -f 1 )
  paste_sample_ids ${opt_hqprefix}.eigenvec \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -u -k 1,1 > ${tmpprefix}.1.bio
cp ${tmpprefix}.1.bio ${opt_biofile}

rm ${tmpprefix}*

