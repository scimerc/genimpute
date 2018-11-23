#!/usr/bin/env bash

trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

extract_sample_ids() {
  if [ ${#@} -le 1 ] ; then
    cat "${1:-/dev/stdin}"
  else
    local -r idfile=$1
    shift
    awk -F $'\t' '{
      if ( NR == FNR ) idarr[$1,$2] = 1
      else if ( ($1,$2) in idarr ) print
    }' "${idfile}" "${*:-/dev/stdin}"
  fi
}

export -f extract_sample_ids

#-------------------------------------------------------------------------------

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

# in: (rs, chr, cm, bp[, a[, m]]) sorted on rs; stdout: same
# get unique chr,cm,bp[,a[,m]] entries [sort|uniq]
# sort on rsnumbers corresponding to unique entries
# get complementary set of entries (duplicated vars) [join -v1]
get_variant_info_for_dup_chr_cm_bp_aa_mm() {
  local -r inputfile="$1"
  sort -k 2 ${inputfile} \
    | uniq -u -f 1 \
    | sort -t ' ' -k 1,1 \
    | join -t ' ' -v1 ${inputfile} -
}

export -f get_variant_info_for_dup_chr_cm_bp_aa_mm

#-------------------------------------------------------------------------------

# in: tab-separated plink bim; out: same
make_variant_names_universal_in_bim_file() {
  local -r fn=$1
  awk -F $'\t' '{
    OFS="\t"
    a[1] = $5; a[2] = $6; asort(a)
    $2 = $1":"$4"_"a[1]"_"a[2]
    if ( $2 in catalog ) {
      catalog[$2]++
      $2 = $2"_"catalog[$2]
    } else catalog[$2] = 1
    print;
  }' ${fn} > ${fn}.tmp
  mv ${fn}.tmp ${fn}
}

export -f make_variant_names_universal_in_bim_file

#-------------------------------------------------------------------------------

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

export -f paste_sample_ids

#-------------------------------------------------------------------------------

tabulate() {
  sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;'
}

export -f tabulate

#-------------------------------------------------------------------------------

