#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"

bcftoolsexec=${BASEDIR}/lib/3rd/bcftools
genimputeexec=${BASEDIR}/lib/progs/genimpute.sh

#-------------------------------------------------------------------------------

# input: unprocessed reference VCF files
# output: processed reference VCF files

usage() {
cat << EOF
USAGE: $( basename $0 ) [\
-c <configuration file> \
-o <output prefix='${opt_outprefix_default}'\
] \
<bcf|vcf file(s)>
 NOTE: $( basename $0 ) uses $( basename ${genimputeexec} )'s quality control functionality.
       <configuration file> is $( basename ${genimputeexec} )'s configuration file.
       (see $( basename ${genimputeexec} ) for more details on available options)
EOF
}

declare -a opt_inputfiles=()
declare opt_outprefix_default=refset
declare opt_outprefix=${opt_outprefix_default}
declare opt_config=''

regexchr='\([Cc]\([Hh]*[Rr]\)*\)*'
regexnum='[0-9xyXYmtMT]\{1,2\}'
regexspecx='\([^[:alnum:]]*\([Nn][Oo][Nn]\)*[Pp][Aa][Rr][12]*\)*'
regexspecxnonpar='chr\([xX]\|23\)\([^[:alnum:]]*[Nn][Oo][Nn][Pp][Aa][Rr]\)*\.bcf\.gz$'
regexspecxpar='chr\([xX][yY]*\|23\|25\)\([^[:alnum:]]*[Pp][Aa][Rr][12]*\)*\.bcf\.gz$'

sedcmd='s/^.*[^[:alnum:]]'${regexchr}'\('${regexnum}${regexspecx}'\)[^[:alnum:]].\+$/chr\3/g'

while getopts "c:o:" opt; do
case "${opt}" in
  c)
    opt_config="${OPTARG}"
    ;;
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

declare  outdir
         outdir="$( dirname "${opt_outprefix}" )"
readonly outdir
declare -r tmpprefix="${opt_outprefix}_tmp"

IFS=$'\n' opt_inputfiles=( $( printf "%s\n" "$@" | sort -V ) ) ; unset IFS

if [ "${#opt_inputfiles[@]}" -eq 0 ] ; then
  usage
  exit 1
else
  ls "${opt_inputfiles[@]}" > /dev/null
fi

echo 'finding common ID set..'
for k in ${!opt_inputfiles[@]} ; do
  if [ $k -eq 0 ] ; then
    "${bcftoolsexec}" query -l "${opt_inputfiles[$k]}" | sort -u > "${tmpprefix}.0.ids"
  else
    "${bcftoolsexec}" query -l "${opt_inputfiles[$k]}" | sort -u | join - "${tmpprefix}.1.ids" > "${tmpprefix}.0.ids"
  fi
  mv "${tmpprefix}.0.ids" "${tmpprefix}.1.ids"
done

if [ -s "${opt_outprefix}_all.ids" ] ; then
  cmp "${tmpprefix}.1.ids" "${opt_outprefix}_all.ids" > /dev/null || { echo 'IDs changed. aborting..'; exit 0; }
else
  mv "${tmpprefix}.1.ids" "${opt_outprefix}_all.ids"
fi

echo 'extracting common ID set..'
for k in ${!opt_inputfiles[@]} ; do
  chrtag=$( sed "${sedcmd}" <<<"${opt_inputfiles[$k]##*/}" )
  [ -s "${opt_outprefix}_${chrtag}.bcf.gz" ] && continue
  "${bcftoolsexec}" view -S "${opt_outprefix}_all.ids" "${opt_inputfiles[$k]}" -Ob > "${tmpprefix}_rfn${k}.bcf.gz"
  mv "${tmpprefix}_rfn${k}.bcf.gz" "${opt_outprefix}_${chrtag}.bcf.gz"
done

echo 'merging all chromosomes..'
chrlist=$( find "${outdir}" -name "${opt_outprefix##*/}"'_chr*.bcf.gz' || true )
[ -s "${opt_outprefix}_all.bcf.gz" ] || ${bcftoolsexec} concat ${chrlist} -Ob > "${tmpprefix}_all.bcf.gz"

echo 'merging non-pseudoautosomal X-chunks..'
xnonparlist=$(
  find "${outdir}" -name "${opt_outprefix##*/}"'_chr*.bcf.gz' \
    | grep ${regexspecxnonpar} | grep -v ${regexspecxpar} || true
)
[ -s "${opt_outprefix}_chr23.bcf.gz" ] || ${bcftoolsexec} concat ${xnonparlist} -Ob > "${tmpprefix}_chr23.bcf.gz"

echo 'merging pseudoautosomal X-chunks..'
xparlist=$(
  find "${outdir}" -name "${opt_outprefix##*/}"'_chr*.bcf.gz' \
    | grep ${regexspecxpar} | grep -v ${regexspecxnonpar} || true
)
[ -s "${opt_outprefix}_chr25.bcf.gz" ] || ${bcftoolsexec} concat ${xparlist} -Ob > "${tmpprefix}_chr25.bcf.gz"

"${bcftoolsexec}" query -f '%CHROM\t%ID\t0\t%POS\t%REF\t%ALT\n' "${tmpprefix}_all.bcf.gz" \
  | tee "${tmpprefix}_all.bim" | cut -f 1,4,5,6 > "${tmpprefix}_all.gpa"

if [ ! -s "${opt_outprefix}_all.unrel.ids" ] ; then
  ${genimputeexec} "${opt_config}" -o "${tmpprefix}" -q "${tmpprefix}_all.bcf.gz"
  mv "${tmpprefix}_"*/.i/qc/e_indqc.ids "${opt_outprefix}_all.unrel.ids"
fi

rename "${tmpprefix}" "${opt_outprefix}" "${tmpprefix}_chr"* "${tmpprefix}_all"*

rm -vrf "${tmpprefix}"*

