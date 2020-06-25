#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"

#-------------------------------------------------------------------------------

# input: single genomic map file for use with Eagle
# output: multiple chromosome-wide genomic map files for separate use with Eagle

usage() {
cat << EOF
USAGE: $( basename $0 ) [\
-o <output prefix='${opt_outprefix_default}'\
] \
<genetic map file>
EOF
}

declare -a opt_inputfiles=()
declare opt_outprefix_default=genomemap
declare opt_outprefix=${opt_outprefix_default}

while getopts "o:" opt; do
case "${opt}" in
  o)
    opt_outprefix="${OPTARG}"
    ;;
  *)
    usage
    exit 1
    ;;
esac
done
shift $((OPTIND-1))

declare -r tmpprefix=${opt_outprefix}_tmp

if [ "${#@}" -eq 0 ] ; then
  usage
  exit 1
else
  opt_inputfile=$1
  ls "${opt_inputfile[@]}" > /dev/null
fi

chrlist=$( zcat -f "${opt_inputfile}" | tail -n +2 | cut -d ' ' -f 1 | sort -nu )

for chr in ${chrlist} ; do
  echo "processing chromosome ${chr}.."
  zcat -f "${opt_inputfile}" | awk -v chr=${chr} 'NR == 1 || $1 == chr' \
    | gzip -c > "${tmpprefix}_${chr}.txt.gz"
  if [ "${chr}" == "23" ] ; then
    zcat -f "${opt_inputfile}" | awk -v chr=${chr} 'NR == 1 || $1 == chr {
      if ( NR == 1 || $2 < 2699520 || $2 > 154931044 ) print
    }' | gzip -c > "${tmpprefix}_25.txt.gz"
  fi
done
rename "${tmpprefix}_" "${opt_outprefix}_chr" "${tmpprefix}_"*.txt.gz

