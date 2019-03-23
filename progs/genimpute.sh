#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

#-------------------------------------------------------------------------------

# define default configuration

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source ${BASEDIR}/progs/cfgmgr.sh

cfgvar_init_from_file ${BASEDIR}/lib/data/genimpute_default.cfg

#---------------------------------------------------------------------------------------------------

# define default options and parse command line

declare opt_dryimpute=0
declare opt_minivarset=0
declare opt_outprefixdefault='genimpute'
declare opt_outprefixbase=${opt_outprefixdefault}
declare opt_samplewhitelist=""

usage() {
cat << EOF

USAGE: $( basename $0 ) [OPTIONS] <bed|bcf|vcf file(s)>

  where <bed|bcf|vcf file(s)> are the genotype files to be merged and qc'd.

NOTE:
  <bed> stands for plink binary set identifiers. $( basename $0 ) expects to find
  <bim> and <fam> files on the same path.

OPTIONS:
  -c <config file>      optional configuration file
  -d                    dry imputation: write scripts but do not run them
  -m                    reduce variant set to the minimal one common to all
  -w <sample file>      optional white list of individuals to restrict qc to
  -o <output prefix>    optional output prefix [default: '${opt_outprefixdefault}']
  -h                    show help message


CONFIGURATION:
  the <config file> may contain custom definitions for a number of variables.
  [see default configuration file 'genimpute_default.cfg' for more information]

EOF
}

printlog() {
  local lvl
  lvl=$1
  readonly lvl
  local lvl_max
  lvl_max="$( cfgvar_get log_lvl )"
  readonly lvl_max
  local logfn
  logfn="$( cfgvar_get logfn )"
  readonly logfn
  # print struff to log
  IFS=$'\n'
  while read line; do
    echo $(date) $lvl "${line}" >> ${logfn}
    # check log-level and print if lower or equal
    if [ "$lvl" -le "$lvl_max" ]; then
      echo "${line}"
    fi
  done
  unset IFS
  return 0
}
export -f printlog

while getopts "c:dmo:w:h" opt; do
case "${opt}" in
  c)
    opt_cfgfile="${OPTARG}"
    ;;
  d)
    opt_dryimpute=1
    ;;
  m)
    opt_minivarset=1
    ;;
  o)
    opt_outprefixbase="${OPTARG}"
    ;;
  w)
    opt_samplewhitelist="${OPTARG}"
    ;;
  h)
    usage
    exit 0
    ;;
  *)
    usage
    exit 1
    ;;
esac
done
shift $((OPTIND-1))

declare -r opt_inputfiles="$(
cat << EOF
$*
EOF
)"

if [ -z "${opt_inputfiles}" ] ; then
  echo -e "\nyou may have neglected to provide any input genotype files."
  usage
  exit 1
else
  ls ${opt_inputfiles} > /dev/null
fi

#---------------------------------------------------------------------------------------------------

# set user configuration, if any, and print all information

if [ ! -z "${opt_cfgfile+x}" ] ; then
  if [ ! -s "${opt_cfgfile}" ] ; then
    echo "configuration file '${opt_cfgfile}' is unusable. aborting..\n"
    exit 1
  else
    cfgvar_update_from_file "${opt_cfgfile}"
  fi
fi
# lock variables
cfgvar_setall_readonly

# append configuration tag to output prefix
opt_outprefixbase=${opt_outprefixbase}_$( cfgvar_show_config | md5sum | head -c 6 )

# export input and output names
export opt_inputfiles
export opt_outprefixbase

# initialize log file
> $( cfgvar_get logfn )

{ 
  echo
  echo -e "================================================================================"
  echo -e "$( basename $0 ) -- $(date)"
  echo -e "================================================================================"
  echo -e "\ngenotype files:\n$( ls -1 ${opt_inputfiles} )\n"
  echo -e "================================================================================"
  echo
} | printlog 1

# print variables
cfgvar_show_config | printlog 1

echo -e "\n================================================================================\n" \
  | printlog 1

#---------------------------------------------------------------------------------------------------

# define executables

plinkexec="${BASEDIR}/lib/3rd/plink --allow-extra-chr"
export plinkexec

#---------------------------------------------------------------------------------------------------

# check dependencies

{
  locale
  echo
  awk --version 2> /dev/null || { echo 'awk is not installed. aborting..'; exit 1; }
  echo
  join --version 2> /dev/null || { echo 'join is not installed. aborting..'; exit 1; }
  echo
  R --version 2> /dev/null || { echo 'R is not installed. aborting..'; exit 1; }
  echo -e "================================================================================"
  echo
} | printlog 1

#---------------------------------------------------------------------------------------------------

# source utility functions

source ${BASEDIR}/progs/qc-tools.sh

#---------------------------------------------------------------------------------------------------
# pre-processing
#---------------------------------------------------------------------------------------------------

# export vars
export opt_dryimpute
export opt_minivarset
export opt_samplewhitelist
export opt_varwhitelist=${opt_outprefixbase}_vwlist.mrk
export opt_refprefix=${opt_outprefixbase}_refset

#---------------------------------------------------------------------------------------------------

echo "==== Whitelist ====" | printlog 0

export opt_outprefix=${opt_outprefixbase}_vwlist

# call wlist
bash ${BASEDIR}/progs/qc-wlist.sh

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

echo "==== Recoding ====" | printlog 0

export opt_outprefix=${opt_outprefixbase}_a_recode

# call recode
bash ${BASEDIR}/progs/qc-recode.sh

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

echo "==== Alignment ====" | printlog 0

export opt_inprefix=${opt_outprefixbase}_a_recode
export opt_outprefix=${opt_outprefixbase}_b_align

# call align
bash ${BASEDIR}/progs/qc-align.sh

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# merge batches (if more than a single one)

echo "==== Merge ====" | printlog 0

export opt_inprefix=${opt_outprefixbase}_b_align
export opt_outprefix=${opt_outprefixbase}_c_proc

# call merge
bash ${BASEDIR}/progs/qc-merge.sh

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# initialize qc iteration counter TODO: implement actual QC loop
declare qciter=0

# biography

if [ ! -f ${opt_outprefixbase}.bio ] ; then
  echo "initializing sample biography file.."
  # initialize sample biography file
  declare cfg_uid
  cfg_uid="$( cfgvar_get uid )"; readonly cfg_uid
  declare opt_outprefix=${opt_outprefixbase}_c_proc
  declare opt_biofile=${opt_outprefixbase}.bio
  cut -f 1-4 ${opt_outprefix}.fam \
    | awk -F $'\t' -v uid=${cfg_uid} -f ${BASEDIR}/lib/awk/idclean.awk \
      --source '
        BEGIN{
          OFS="\t";
          print( uid, "FID", "IID", "P1ID", "P2ID" )
        } { print( idclean( $1"_"$2 ), $0 ) } \
      ' \
    | sort -u -k 1,1 > ${opt_biofile}
  Norg=$( cat ${opt_outprefix}.fam | wc -l )
  Nuid=$( cat ${opt_biofile} | wc -l )
  if [ ${Norg} -ge ${Nuid} ] ; then
    echo '========> conflicts generated in universal IDs.'
    echo 'please recode IDs so they do not lead to conflicts.'
    exit 1
  fi
  unset opt_outprefix
  unset opt_biofile
fi

#---------------------------------------------------------------------------------------------------

# get high quality set

echo "==== Variant HQC ====" | printlog 0

# export vars
export opt_inprefix=${opt_outprefixbase}_c_proc
export opt_hqprefix=${opt_outprefixbase}_d_hqset

# call hqset
bash ${BASEDIR}/progs/qc-hqset.sh

# cleanup
unset opt_hqprefix
unset opt_inprefix

#---------------------------------------------------------------------------------------------------

# identify duplicate, mixup and related individuals

echo "==== Individual QC ====" | printlog 0

# export vars
export opt_inprefix=${opt_outprefixbase}_c_proc
export opt_hqprefix=${opt_outprefixbase}_d_hqset
export opt_outprefix=${opt_outprefixbase}_e_indqc
export opt_biofile=${opt_outprefixbase}.bio

# call indqc
bash ${BASEDIR}/progs/qc-ind.sh

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

echo "==== Variant QC ====" | printlog 0

# export vars
export opt_hqprefix=${opt_outprefixbase}_d_hqset
export opt_inprefix=${opt_outprefixbase}_e_indqc
export opt_outprefix=${opt_outprefixbase}_f_varqc

# call qcvar
bash ${BASEDIR}/progs/qc-var.sh

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# perform final QC (control-HWE? batch effects?)

# export vars
export opt_hqprefix=${opt_outprefixbase}_d_hqset
export opt_inprefix=${opt_outprefixbase}_f_varqc
export opt_outprefix=${opt_outprefixbase}_g_finqc
export opt_batchoutprefix=${opt_outprefixbase}_c_proc_batch

# call qcfinal
bash ${BASEDIR}/progs/qc-final.sh

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_batchoutprefix

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix=${opt_outprefixbase}_g_finqc
export opt_hqprefix=${opt_outprefixbase}_h_hqset

# call hqset
bash ${BASEDIR}/progs/qc-hqset.sh $( cfgvar_get refprefix )

# cleanup
unset opt_hqprefix
unset opt_inprefix

#---------------------------------------------------------------------------------------------------

# compute genetic PCs

echo "==== PCA etc. ====" | printlog 0

# export vars
export opt_hqprefix=${opt_outprefixbase}_h_hqset
export opt_inprefix=${opt_outprefixbase}_g_finqc
export opt_outprefix=${opt_outprefixbase}_g_finqc
export opt_biofile=${opt_outprefixbase}.bio

# call getpcs
bash ${BASEDIR}/progs/qc-getpcs.sh

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# convert

echo "==== Recoding ====" | printlog 0

# export vars
export opt_inprefix=${opt_outprefixbase}_g_finqc
export opt_outprefix=${opt_outprefixbase}_i_pre

# call convert
bash ${BASEDIR}/progs/convert.sh

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# impute

# export vars
export opt_inprefix=${opt_outprefixbase}_i_pre
export opt_outprefix=${opt_outprefixbase}_i_imp

# call impute
bash ${BASEDIR}/progs/impute.sh

# cleanup
unset opt_inprefix
unset opt_outprefix

