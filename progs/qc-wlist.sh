#!/usr/bin/env bash

# exit on error
trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -r  tmpprefix=${opt_outprefix}_tmp
declare -r  debuglogfn=${tmpprefix}_debug.log
declare -ra batchfiles=( ${opt_inputfiles} )

declare -r cfg_genomebuild="$( cfgvar_get genomebuild )"
declare -r cfg_refallelesfn="$( cfgvar_get refallelesfn )"


if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping final QC step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
  exit 1
fi

for i in ${!batchfiles[@]} ; do
  declare plinkflag=""
  case "$( get_genotype_file_format "${batchfiles[$i]}" )" in
    "bed" ) 
      plinkflag="--bim"
      ;;
    "vcf" )
      plinkflag="--vcf"
      ;;
    "bcf" )
      plinkflag="--bcf"
      ;;
    * ) 
      printf "error: fileformat not handled\n" >&2
      exit 1
      ;;
  esac
  # whatever format the input file is - make a bim file
  plink ${plinkflag} ${batchfiles[$i]/%.bed/.bim} \
        --make-just-bim \
        --out ${tmpprefix}_ex \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  sed -i -r 's/[ \t]+/\t/g; s/^chr//g;' ${tmpprefix}_ex.bim
  sed -i -r 's/^XY/X/g; s/^X/23/g; s/^Y/24/g; s/^25/23/g;' ${tmpprefix}_ex.bim
  if [ -s "${tmpprefix}_ex.bim" ] ; then
    if [ $i -eq 0 ] ; then
      bimtogprs ${tmpprefix}_ex.bim \
        | sort -t $'\t' -u -k 1,1 \
        > ${tmpprefix}_ex.0.gprs
    else
      join -t $'\t' ${tmpprefix}_ex.1.gprs <( bimtogprs ${tmpprefix}_ex.bim ) \
        | cut -f 1,2 \
        | sort -t $'\t' -u -k 1,1 \
        > ${tmpprefix}_ex.0.gprs
    fi
    mv ${tmpprefix}_ex.0.gprs ${tmpprefix}_ex.1.gprs
    rm -f ${tmpprefix}_ex.bim
  fi
  awk -F $'\t' '{
    OFS="\t"
    split( $1, gpvec, ":" )
    print( gpvec[1], gpvec[2] - 1, gpvec[2], $2 )
  }' ${tmpprefix}_ex.1.gprs > ${opt_varwhitelist}
  unset plinkflag
done


