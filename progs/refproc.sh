#!/usr/bin/env bash
#TODO: generate m3vcf files

# exit on error
set -ETeuo pipefail

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"

# uncomment this for dry run
# declare -r DRYMODE='drymode'

bcftoolsexec=${BASEDIR}/lib/3rd/bcftools
genimputeexec=${BASEDIR}/progs/genimpute.sh

#-------------------------------------------------------------------------------

# input: unprocessed reference VCF files
# output: processed reference VCF files

usage() {
cat << EOF

$( basename $0 ) pre-processes reference data for use in $( basename ${genimputeexec} ).

please note that the Y chromosome is currently excluded from parts of the processing.
imputation of the Y chromosome will have to be performed separately.

USAGE:
  $( basename $0 ) [\
-c <configuration file> \
-o <output prefix='${opt_outprefix_default}>'\
] \
<bcf|vcf file(s)>

 NOTE:
  $( basename $0 ) uses $( basename ${genimputeexec} )'s own quality control functionality.
  <configuration file> is $( basename ${genimputeexec} )'s configuration file.
  (see $( basename ${genimputeexec} ) for more details on available options)

EOF
}

declare -a opt_inputfiles=()
declare opt_outprefix_default=refset
declare opt_outprefix=${opt_outprefix_default}
declare opt_config=''

Nmin=5
regexchr='\([Cc]\([Hh]*[Rr]\)*\)*'
regexnum='[0-9xyXYmtMT]\{1,2\}'
regexspecx='\([^[:alnum:]]*\([Nn][Oo][Nn]\)*[^[:alnum:]]*[Pp][AaSs][RrEe][Uu]*[Dd]*[Oo]*[12]*\)*'
regexspecxnonpar='chr\([xX]\|23\)\([^[:alnum:]]*[Nn][Oo][Nn][^[:alnum:]]*[Pp][AaSs][RrEe][Uu]*[Dd]*[Oo]*\)*\.bcf\.gz$'
regexspecxpar='chr\([xX][yY]*\|23\|25\)\([^[:alnum:]]*[Pp][AaSs][RrEe][Uu]*[Dd]*[Oo]*[12]*\)*\.bcf\.gz$'
regexspecy='chr\([yY]\|24\)\.bcf\.gz$'

grepcmd=${regexchr}'\('${regexnum}${regexspecx}'\)'

gpaquery='%CHROM\t%ID\t0\t%POS\t%REF\t%ALT\n'

while getopts "c:o:" opt; do
case "${opt}" in
  c)
    opt_config="-c ${OPTARG}"
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

dflag=0
echo 'finding common ID set..'
for k in ${!opt_inputfiles[@]} ; do
  if [ $k -eq 0 ] ; then
    "${bcftoolsexec}" query -l "${opt_inputfiles[$k]}" | sort -u > "${tmpprefix}.0.ids"
  else
    "${bcftoolsexec}" query -l "${opt_inputfiles[$k]}" | sort -u | join - "${tmpprefix}.1.ids" > "${tmpprefix}.0.ids"
    cmp "${tmpprefix}.0.ids" "${tmpprefix}.1.ids" > /dev/null || dflag=1
  fi
  mv "${tmpprefix}.0.ids" "${tmpprefix}.1.ids"
done

# compare to existing common ID set if one is there
if [ -s "${opt_outprefix}_all.ids" ] ; then
  cmp "${tmpprefix}.1.ids" "${opt_outprefix}_all.ids" > /dev/null || { echo 'IDs changed. aborting..'; exit 0; }
else
  mv "${tmpprefix}.1.ids" "${opt_outprefix}_all.ids"
fi

echo 'processing files..'
if [ "${DRYMODE:-}" == "" ] ; then
  for k in ${!opt_inputfiles[@]} ; do
    chrtag=$( grep -o "${grepcmd}" <<<"${opt_inputfiles[$k]##*/}" | head -n 1 || true )
    [ -s "${opt_outprefix}_${chrtag}.bcf.gz" ] && continue
    # extract common ID set if any differences were detected
    echo "processing chromosome '${chrtag}'.."
    samplefilter=''
    [ ${dflag} -eq 1 ] && samplefilter="-S ${opt_outprefix}_all.ids" 
    "${bcftoolsexec}" view ${samplefilter} -c ${Nmin}:minor "${opt_inputfiles[$k]}" -Ob \
      > "${tmpprefix}_rfn${k}.bcf.gz"
    "${bcftoolsexec}" index "${tmpprefix}_rfn${k}.bcf.gz"
    rename "${tmpprefix}_rfn${k}" "${opt_outprefix}_${chrtag}" "${tmpprefix}_rfn${k}".*
  done
fi

bcftoolscmd='convert'
chrlist=( $(
  find "${outdir}" -maxdepth 1 -name "${opt_outprefix##*/}"'_chr*.bcf.gz' \
    | grep -v ${regexspecy}
) )
[ ${#chrlist[@]} -gt 1 ] && bcftoolscmd='concat -a'
cmdall="${bcftoolsexec} ${bcftoolscmd} ${chrlist[@]} -Ob"
bcftoolscmd='convert'
xnonparlist=( $(
  find "${outdir}" -maxdepth 1 -name "${opt_outprefix##*/}"'_chr*.bcf.gz' \
    | grep ${regexspecxnonpar} | grep -v -e ${regexspecxpar} -e ${regexspecy} || true
) )
[ ${#xnonparlist[@]} -gt 1 ] && bcftoolscmd='concat -a'
cmdxnonpar="${bcftoolsexec} ${bcftoolscmd} ${xnonparlist[@]} -Ob"
bcftoolscmd='convert'
xparlist=( $(
  find "${outdir}" -maxdepth 1 -name "${opt_outprefix##*/}"'_chr*.bcf.gz' \
    | grep ${regexspecxpar} | grep -v -e ${regexspecxnonpar} -e ${regexspecy} || true
) )
[ ${#xparlist[@]} -gt 1 ] && bcftoolscmd='concat -a'
cmdxpar="${bcftoolsexec} ${bcftoolscmd} ${xparlist[@]} -Ob"
bcftoolscmd='convert'
ylist=( $(
  find "${outdir}" -maxdepth 1 -name "${opt_outprefix##*/}"'_chr*.bcf.gz' \
    | grep ${regexspecy} || true
) )
[ ${#ylist[@]} -gt 1 ] && bcftoolscmd='concat -a'
cmdy="${bcftoolsexec} ${bcftoolscmd} ${ylist[@]} -Ob"
if [ "${DRYMODE:-}" == "" ] ; then
  if [ ! -s "${opt_outprefix}_all.bcf.gz" ] ; then
    echo 'merging all chromosomes (except Y)..'
    $cmdall > "${tmpprefix}_all.bcf.gz"
    "${bcftoolsexec}" index "${tmpprefix}_all.bcf.gz"
    rename "${tmpprefix}" "${opt_outprefix}" "${tmpprefix}_all.bcf.gz"*
  fi
  if [ ! -s "${opt_outprefix}_chr23.bcf.gz" ] ; then
    echo 'merging non-par X regions..'
    $cmdxnonpar > "${tmpprefix}_chr23.bcf.gz"
    "${bcftoolsexec}" index "${tmpprefix}_chr23.bcf.gz"
    rename "${tmpprefix}" "${opt_outprefix}" "${tmpprefix}_chr23.bcf.gz"*
  fi
  if [ ! -s "${opt_outprefix}_chr25.bcf.gz" ] ; then
    echo 'merging par X regions..'
    $cmdxpar > "${tmpprefix}_chr25.bcf.gz"
    "${bcftoolsexec}" index "${tmpprefix}_chr25.bcf.gz"
    rename "${tmpprefix}" "${opt_outprefix}" "${tmpprefix}_chr25.bcf.gz"*
  fi
  if [ ! -s "${opt_outprefix}_chr24.bcf.gz" ] ; then
    echo 'merging par Y regions..'
    $cmdy > "${tmpprefix}_chr24.bcf.gz"
    "${bcftoolsexec}" index "${tmpprefix}_chr24.bcf.gz"
    rename "${tmpprefix}" "${opt_outprefix}" "${tmpprefix}_chr24.bcf.gz"*
  fi
  if [ ! -s "${opt_outprefix}_all.gpa" ] ; then
    echo 'generating variant list..'
    "${bcftoolsexec}" query -f "${gpaquery}" "${opt_outprefix}_all.bcf.gz" \
      | tee "${tmpprefix}_all.bim" | cut -f 1,4,5,6 > "${tmpprefix}_all.gpa"
    rename "${tmpprefix}" "${opt_outprefix}" "${tmpprefix}_all"*
  fi
  if [ ! -s "${opt_outprefix}_all.unrel.ids" ] ; then
    echo 'finding unrelated individuals..'
    ${genimputeexec} ${opt_config} -q -o "${tmpprefix}" "${opt_outprefix}_all.bcf.gz"
    mv "${tmpprefix}_"*/.i/qc/e_indqc.ids "${opt_outprefix}_all.unrel.ids"
  fi
else
  echo $cmdall
  echo $cmdxnonpar
  echo $cmdxpar
  echo $cmdy
fi

rm -vrf "${tmpprefix}"*

