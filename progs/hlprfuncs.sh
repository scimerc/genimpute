#!/usr/bin/env bash

function bdy() {
  cnt=1
  unset myargs
  mycmd=$1; shift
  for carg in $* ; do
    if [[ -f "${carg}" ]] ; then
      myargs[${cnt}]="${carg}"
      cnt=$((cnt + 1)) 
    else
      mycmd="${mycmd} ${carg}"
    fi
  done
  cnt=1
  if (( ${#myargs[*]} > 0 )) ; then
    if (( cnt == 1 )) ; then
      head -n 1 ${myargs[${cnt}]}
      tail -n +2 ${myargs[${cnt}]} | ${mycmd}
    else
      ${mycmd} ${myargs[${cnt}]}
    fi
    cnt=$((cnt + 1)) 
  else
    read myheader
    printf "%s\n" "${myheader}"
    ${mycmd} /dev/stdin
  fi
}

export -f bdy

extract_sample_ids() {
  local idfile
  idfile="$1"
  readonly idfile
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
  local -r inputfile="$1"
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
    printf "> error: get_genotype_file_format - file '%s' not found\n" "${inputfile}" >&2
    return 1
  fi

  local    filehead
           filehead=$( zcat -f "${inputfile}" | hexdump -n2 -e '2/1 "0x%02x "' )
  readonly filehead
  case ${filehead} in
    "${bedhex}" )
      printf "bed" ;;
    "${bcfhex}" )
      printf "bcf" ;;
    "${vcfhex}" )
      printf "vcf" ;;
    * )
      printf "> error: unknown file format '%s'\n" "${filehead}" >&2
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
           md5=$(echo "${path}" | md5sum)
  readonly md5
  printf "%s_%s" "$(basename "${path%.bed}" .gz)" "${md5:0:$numchars}"
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
  awk -F $'\t' '$1 == 23' "${@:-/dev/stdin}" | wc -l
}

export -f get_xvar_count

#-------------------------------------------------------------------------------

# in: plink bim; out: same
standardize_bim_file() {
  local -r fn="$1"
  awk '{
    OFS="\t"
    a[1] = toupper($5)
    a[2] = toupper($6)
    if ( a[1] <= 0 ) a[1] = 0
    if ( a[2] <= 0 ) a[2] = 0
    gsub( "^[Cc]([Hh][Rr]*)*", "", $1 )
    if ( $1 ~ "[Xx][Yy]" ) $1 == 25
    if ( $1 ~ "^[Xx]"    ) $1 == 23
    if ( $1 ~ "^[Yy]"    ) $1 == 24
    if ( $1 ~ "^[Mm]"    ) $1 == 26
    if ( $1 != 0 && $4 != 0 )
      $2 = $1":"$4":"a[1]":"a[2]
    if ( $2 in catalog ) {
      catalog[$2]++
      $2 = $2"_"catalog[$2]
    } else catalog[$2] = 1
    print;
  }' "${fn}" > "${fn}.tmp"
  mv "${fn}.tmp" "${fn}"
}

export -f standardize_bim_file

#-------------------------------------------------------------------------------

# in: headed file with leading <fid iid> fields (e.g. plink fam file);
# stdout: tab-separated headed file with leading <uid=fid_iid> field, sorted on the latter;
attach_uids() {
  local cmdarg=''
  local header=0
  local infile=''
  local locuid='000UID'
  cmdarg=$1
  readonly cmdarg
  if [ "${cmdarg}" == "-h" ] ; then
    header=1
    readonly header
    shift
  fi
  infile="$1"
  readonly infile
  if [ -s "${infile}" ] ; then
    tabulate "${infile}" \
      | awk -F $'\t' \
        -f "${BASEDIR}/lib/awk/idclean.awk" \
        -v head=${header} \
        -v uid=${locuid} \
        --source '{
          # if the file has no header begin processing from the first line
          if ( head == 0 || NR > 1 )
            uid=idclean($1"_"$2)
          printf( "%s", uid )
          for ( k=1; k<=NF; k++ )
            printf( "\t%s", $k )
          printf( "\n" )
        }'
  fi
}

export -f attach_uids

#-------------------------------------------------------------------------------

tabulate() {
  sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;' "${@:-/dev/stdin}"
}

export -f tabulate

#-------------------------------------------------------------------------------

rename() {
  local -r fromtext=$1
  local -r totext=$2
  shift 2
  for f in "$@" ; do
    [ -e "$f" ] && mv $f ${f/${fromtext}/${totext}}
  done
}

export -f rename

#-------------------------------------------------------------------------------

# using () instead of {} for function body to limit scope of inner function
split_sex() (
  local -r fn_updatesex=$1
  local -r max_groupsize=$2
  local -r outprefix=$3
  local -r tmppath=$4
  local -r sqlite3_bin="${BASEDIR}"/lib/3rd/sqlite3
  local -r sqlite3_db="${tmppath}"_genimpute.db
  
  # local function
  _get_unique_subjids_by_sexid() {
    local -r sqlite3_bin=$1
    local -r sqlite3_db=$2
    local -r sex_id=$3

  sql=$( cat << EOI
.mode "csv"
.separator "\t"
select distinct(samples.subjid)
from samples
left join sex_info
on sex_info.subjid = samples.subjid 
where 
  case sex_info.sex 
    when 1 then 1 
    when 2 then 2 
    else 0         -- NOTE: '0' is unknown sex
  end = ${sex_id}
;
EOI
  )
    ${sqlite3_bin} --init <(echo "${sql}") "${sqlite3_db}"
  }

  # read input table (all samples) from stdin 
  sql=$( cat << EOI
drop table if exists samples;
create table samples (
  subjid text
);
.mode "csv"
.separator " "
.import /dev/stdin samples
EOI
  )
  ${sqlite3_bin} --init <(echo "${sql}") "${sqlite3_db}"
  
  # read sex info; 3 columns, space-separated 
  # example of line: "HG00097_HG00097 HG00097_HG00097 0")
  sql=$( cat << EOI
drop table if exists sex_info;
create table sex_info (
  subjid text
  , fam_id text  --NOTE: this second ID is supposed to be identical to first ID
  , sex int
);
.mode "csv"
.separator " "
.import ${fn_updatesex} sex_info
EOI
  )
  ${sqlite3_bin} --init <(echo "${sql}") "${sqlite3_db}"
  
  # append sex ID to fn_outprefix; ´split´ appends more digits
  local -r sex_ids="0 1 2"
  for i in ${sex_ids}; do
    _get_unique_subjids_by_sexid ${sqlite3_bin} ${sqlite3_db} $i | split -d -l ${max_groupsize} /dev/stdin "${outprefix}"$i
  done

  rm -rf ${sqlite3_db}
)

export -f split_sex
