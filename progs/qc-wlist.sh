#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare  cfg_refprefix
         cfg_refprefix="$( cfgvar_get refprefix )"
readonly cfg_refprefix

if [ ! -z "${cfg_refprefix}" ] ; then
  declare -a batchfiles=( ${opt_inputfiles} "${cfg_refprefix}.all.haplotypes.bcf.gz" )
else
  declare -a batchfiles=( ${opt_inputfiles} )
fi
readonly batchfiles

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

#-------------------------------------------------------------------------------

# input: merged plink set
# output: hq plink set (with imputed sex, if possible)

printf "  * Compile whitelist of variants common to all batches (if enabled)\n" | printlog 1

if [ ${opt_minivarset} -eq 1 ] ; then
  exit 0
fi

if [ -f "${opt_varwhitelist}" ] ; then
  printf "variant white list found. skipping compilation..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

#-------------------------------------------------------------------------------

bimtogprs() {
# in: tab-separated plink bim; stdout: chr:bp, rs
  local -r inputfile="$1"
  awk -F $'\t' '{
    OFS="\t"
    print( $1":"$4, $2 )
  }' ${inputfile} \
  | sort -u -k 1,1
}

#-------------------------------------------------------------------------------

for i in ${!batchfiles[@]} ; do
  declare plinkflag=""
  [ -f "${batchfiles[$i]}" ] || { printf "file '%s' not found." "${batchfiles[$i]}"; exit 1; }
  declare fformat=$( get_genotype_file_format "${batchfiles[$i]}" )
  case "${fformat}" in
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
      printf "error: fileformat '%s' not handled\n" ${fformat} >&2
      exit 1
      ;;
  esac
  # whatever format the input file is - make a bim file
  printf "* Recode variants set into bim format" | printlog 1
  ${plinkexec} ${plinkflag} ${batchfiles[$i]/%.bed/.bim} \
        --make-just-bim \
        --out ${tmpprefix}_ex \
        2>&1 | printlog 2
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

rm ${tmpprefix}*

