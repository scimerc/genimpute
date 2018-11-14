#!/usr/bin/env bash

trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

get_genotype_file_format() {
  local -r inputfile=$1
  #TODO implement plink2 native .gen format
  local -r bedhex='0x6c 0x1b'
  local -r bcfhex='0x42 0x43'
  local -r vcfhex='0x23 0x23'

  # get first 2 bytes of inputfile
  if [ -z "${inputfile}" -o ! -f "${inputfile}" ]; then
    return 1
  fi

  local -r filehead=$( zcat -f ${inputfile} | hexdump -n2 -e '2/1 "0x%02x "' )
  case ${filehead} in
    "${bedhex}" )
      printf 'bed' ;;
    "${bcfhex}" )
      printf 'bcf' ;;
    "${vcfhex}" )
      printf 'vcf' ;;
    * )
      printf "error: unknown file format '%s'\n" "${filehead}" >&2
      return 1 ;;
  esac
  return 0
}

export -f get_genotype_file_format

#-------------------------------------------------------------------------------

get_unique_filename_from_path() {
  local -r path="$1"
  local -r numchars=6
  local -r md5=$(echo "$path" | md5sum)
  printf "%s_%s" "$(basename "$path" )" "${md5:0:$numchars}"
}

export -f get_unique_filename_from_path

#-------------------------------------------------------------------------------

# in: tab-separated plink bim; out: same
make_variant_names_universal_in_bim_file() {
  local -r fn=$1
  awk -F $'\t' '{
    OFS="\t";
    a[1] = $5; a[2] = $6; asort(a);
    $2 = $1":"$4"_"a[1]"_"a[2];
    print;
  }' ${fn} > ${fn}.tmp
  mv ${fn}.tmp ${fn}
}

export -f make_variant_names_universal_in_bim_file

#-------------------------------------------------------------------------------

