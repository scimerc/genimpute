#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

#---------------------------------------------------------------------------------------------------

# define default configuration

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source "${BASEDIR}/progs/cfgmgr.sh"

cfgvar_init_from_file "${BASEDIR}/lib/genimpute_default.cfg"

#---------------------------------------------------------------------------------------------------

# define default options and parse command line

declare opt_dryimpute=0
declare opt_minivarset=0
declare opt_qconly=0
declare opt_outprefixdefault='genimpute'
declare opt_outprefixbase="${opt_outprefixdefault}"
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
  -q                    perform quality control only, do not phase or impute
  -w <sample file>      optional white list of individuals to restrict qc to
  -o <output prefix>    optional output prefix [default: '${opt_outprefixdefault}']
  -h                    show help message


CONFIGURATION:
  the <config file> may contain custom definitions for a number of variables.
  [see default configuration file 'genimpute_default.cfg' for more information]

EOF
}

while getopts "c:dmqw:o:h" opt; do
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
  q)
    opt_qconly=1
    ;;
  w)
    opt_samplewhitelist="${OPTARG}"
    ;;
  o)
    opt_outprefixbase="${OPTARG}"
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
opt_outprefixbase="${opt_outprefixbase}_$( cfgvar_show_config | md5sum | head -c 6 )"

# prepare output containers
mkdir -p "${opt_outprefixbase}/bcf" "${opt_outprefixbase}/bed" "${opt_outprefixbase}/.i/.s" \
         "${opt_outprefixbase}/.i/qc" "${opt_outprefixbase}/.i/ph" "${opt_outprefixbase}/.i/im"

# export input and output names
export opt_inputfiles
export opt_outprefixbase

printlog() {
  local lvl
  lvl=$1
  readonly lvl
  local lvl_max
  lvl_max="$( cfgvar_get log_lvl )"
  readonly lvl_max
  local logfn
  logfn="${opt_outprefixbase}.log"
  readonly logfn
  # print stuff to log
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

# initialize log file
> "${opt_outprefixbase}.log"

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

source "${BASEDIR}/progs/qc-tools.sh"

#---------------------------------------------------------------------------------------------------
# pre-processing
#---------------------------------------------------------------------------------------------------

echo "starting.."

# export vars
export cfg_uid
cfg_uid="$( cfgvar_get uid )"
readonly cfg_uid
export cfg_uid
export opt_dryimpute
export opt_minivarset
export opt_samplewhitelist
export opt_varwhitelist="${opt_outprefixbase}/.i/qc/a_vwlist.mrk"
export opt_refprefix="${opt_outprefixbase}/.i/qc/a_refset"

#---------------------------------------------------------------------------------------------------

echo "==== Whitelist ====" | printlog 1

export opt_outprefix="${opt_outprefixbase}/.i/qc/a_vwlist"

# call wlist
bash "${BASEDIR}/progs/qc-wlist.sh"

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

echo "==== Recoding ====" | printlog 1

export opt_outprefix="${opt_outprefixbase}/.i/qc/a_recode"

# call recode
bash "${BASEDIR}/progs/qc-recode.sh"

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

echo "==== Alignment ====" | printlog 1

export opt_inprefix="${opt_outprefixbase}/.i/qc/a_recode"
export opt_outprefix="${opt_outprefixbase}/.i/qc/b_align"

# call align
bash "${BASEDIR}/progs/qc-align.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# merge batches (if more than a single one)

echo "==== Merge ====" | printlog 1

export opt_inprefix="${opt_outprefixbase}/.i/qc/b_align"
export opt_outprefix="${opt_outprefixbase}/.i/qc/c_proc"

# call merge
bash "${BASEDIR}/progs/qc-merge.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# initialize qc iteration counter TODO: implement actual QC loop
declare qciter=0

# biography

# export vars
export opt_outprefix="${opt_outprefixbase}/.i/qc/c_proc"
export opt_biofile="${opt_outprefixbase}/bio.txt"

# call init-bio
bash "${BASEDIR}/progs/init-bio.sh"

# cleanup
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# get high quality set

echo "==== Variant HQC ====" | printlog 1

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/c_proc"
export opt_hqprefix="${opt_outprefixbase}/.i/qc/d_hqset"

# call hqset
bash "${BASEDIR}/progs/qc-hqset.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix

#---------------------------------------------------------------------------------------------------

# identify duplicate, mixup and related individuals

echo "==== Individual QC ====" | printlog 1

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/c_proc"
export opt_hqprefix="${opt_outprefixbase}/.i/qc/d_hqset"
export opt_outprefix="${opt_outprefixbase}/.i/qc/e_indqc"
export opt_biofile="${opt_outprefixbase}/bio.txt"

# call indqc
bash "${BASEDIR}/progs/qc-ind.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

echo "==== Variant QC ====" | printlog 1

# export vars
export opt_hqprefix="${opt_outprefixbase}/.i/qc/d_hqset"
export opt_inprefix="${opt_outprefixbase}/.i/qc/e_indqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/f_varqc"

# call qcvar
bash "${BASEDIR}/progs/qc-var.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# perform final QC (control-HWE? batch effects?)

# export vars
export opt_hqprefix="${opt_outprefixbase}/.i/qc/d_hqset"
export opt_inprefix="${opt_outprefixbase}/.i/qc/f_varqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_batchoutprefix="${opt_outprefixbase}/.i/qc/c_proc_batch"

# call qcfinal
bash "${BASEDIR}/progs/qc-final.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_batchoutprefix

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_hqprefix="${opt_outprefixbase}/.i/qc/h_hqset"

# call hqset
bash "${BASEDIR}/progs/qc-hqset.sh" "$( cfgvar_get refprefix )"

# cleanup
unset opt_hqprefix
unset opt_inprefix

#---------------------------------------------------------------------------------------------------

# compute genetic PCs

echo "==== PCA etc. ====" | printlog 1

# export vars
export opt_hqprefix="${opt_outprefixbase}/.i/qc/h_hqset"
export opt_inprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_biofile="${opt_outprefixbase}/bio.txt"

# call getpcs
bash "${BASEDIR}/progs/qc-getpcs.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

if [ ${opt_qconly} -eq 1 ] ; then
  echo "done."
  exit 0
fi

# convert

echo "==== Recoding ====" | printlog 1

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_outprefix="${opt_outprefixbase}/.i/ph/i_pre"

# call convert
bash "${BASEDIR}/progs/convert.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# impute

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/ph/i_pre"
export opt_outprefix="${opt_outprefixbase}/.i/im/i_imp"

# call impute
bash "${BASEDIR}/progs/impute.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

echo "done."

