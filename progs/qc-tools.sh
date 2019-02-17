#!/usr/bin/env bash

extract_sample_ids() {
  local -r idfile=$1
  if [ -z "${idfile}" ] ; then
    cat
  else
    awk '{
      if ( NR == FNR ) idarr[$1,$2] = 1
      else if ( ($1,$2) in idarr ) print
    }' "${idfile}" /dev/stdin
  fi
}

export -f extract_sample_ids

#-------------------------------------------------------------------------------

get_genotype_file_format() {
  local -r inputfile=$1
  #TODO implement plink2 native .gen format
  local    bedhex
           bedhex='0x6c 0x1b'
  readonly bedhex        
  local    bcfhex
           bcfhex='0x42 0x43'
  readonly bcfhex
  local    vcfhex
           vcfhex='0x23 0x23'
  readonly vcfhex

  # get first 2 bytes of inputfile
  if [ -z "${inputfile}" -o ! -f "${inputfile}" ]; then
    printf "error: get_genotype_file_format - file '%s' not found\n" "${inputfile}" >&2
    return 1
  fi

  local    filehead
           filehead=$( zcat -f ${inputfile} | hexdump -n2 -e '2/1 "0x%02x "' )
  readonly filehead
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
  local    md5
           md5=$(echo "$path" | md5sum)
  readonly md5
  printf "%s_%s" "$(basename "$path" )" "${md5:0:$numchars}"
}

export -f get_unique_filename_from_path

#-------------------------------------------------------------------------------

# in: space-separated (rs, chr, cm, bp[, a[, m]]) sorted on rs; stdout: same
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

# in: plink bim
get_xvar_count() {
  awk '$1 == 23' $1 | wc -l
}

export -f get_xvar_count

#-------------------------------------------------------------------------------

# in: tab-separated plink bim; out: same
standardize_bim_file() {
  local -r fn=$1
  awk -F $'\t' '{
    OFS="\t"
    a[1] = toupper($5)
    a[2] = toupper($6)
    if ( a[1] <= 0 ) a[1] = 0
    if ( a[2] <= 0 ) a[2] = 0
    asort(a)
    if ( $1 ~ "^X" )  $1 == 23
    if ( $1 ~ "^XY" ) $1 == 25
    if ( $1 ~ "^Y" )  $1 == 24
    if ( $1 ~ "^MT" ) $1 == 26
    $2 = $1":"$4"_"a[1]"_"a[2]
    if ( $2 in catalog ) {
      catalog[$2]++
      $2 = $2"_"catalog[$2]
    } else catalog[$2] = 1
    print;
  }' ${fn} > ${fn}.tmp
  mv ${fn}.tmp ${fn}
}

export -f standardize_bim_file

#-------------------------------------------------------------------------------

paste_sample_ids() {
  local -r infile="$1"
  if [ -s "${infile}" ] ; then
    tabulate "${infile}" \
      | awk -F $'\t' -v uid=${cfg_uid} '{
          OFS="\t"
          if ( NR>1 ) uid=$1"_"$2
          printf( "%s", uid )
          for ( k=3; k<=NF; k++ )
            printf( "\t%s", $k )
          printf( "\n" )
        }' \
      | sort -t $'\t' -u -k 1,1
  fi
}

export -f paste_sample_ids

#-------------------------------------------------------------------------------

tabulate() {
  sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;' "${*:-/dev/stdin}"
}

export -f tabulate

#-------------------------------------------------------------------------------
