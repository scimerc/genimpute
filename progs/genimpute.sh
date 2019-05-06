#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

#---------------------------------------------------------------------------------------------------

echo "starting.."

# define default configuration

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source "${BASEDIR}/progs/cfgmgr.sh"

cfgvar_init_from_file "${BASEDIR}/lib/genimpute_default.cfg"

#---------------------------------------------------------------------------------------------------

# define default options and parse command line

declare opt_dryimpute=0
declare opt_imputeonly=false
declare opt_qconly=false
declare opt_outprefixdefault='genimpute'
declare opt_outprefixbase="${opt_outprefixdefault}"

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
  -i                    perform imputation only (assumes qc'd files present)
  -q                    perform quality control only, do not phase or impute
  -o <output prefix>    optional output prefix [default: '${opt_outprefixdefault}']
  -h                    show this help message and exit


CONFIGURATION:
  the <config file> may contain custom definitions for a number of variables.
  [see default configuration file 'genimpute_default.cfg' for more information]

EOF
}

while getopts "c:diqo:h" opt; do
case "${opt}" in
  c)
    opt_cfgfile="${OPTARG}"
    ;;
  d)
    opt_dryimpute=1
    ;;
  i)
    opt_imputeonly=true
    ;;
  q)
    opt_qconly=true
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
opt_outprefixbase="${opt_outprefixbase}_$(
  cfgvar_show_config | egrep -v 'execommand|plink' \
    | md5sum | head -c6
)"

# export input and output names
export opt_inputfiles
export opt_outprefixbase

# prepare output containers
mkdir -p "${opt_outprefixbase}/bcf" "${opt_outprefixbase}/bed" \
         "${opt_outprefixbase}/.i/.s" "${opt_outprefixbase}/.i/qc" "${opt_outprefixbase}/.i/im"

#---------------------------------------------------------------------------------------------------

# define executables

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

declare jobdepflag=""
declare joblist=""

#TODO: find the time executable in the cluster and restore the time commands
# declare -r timeformat="ResStats\tRealSec:\t%e\tMaxMemKB:\t%M\tUserSec:\t%U\tSysSec:\t%S"
# declare -r timexec="/usr/bin/env time -o /dev/stdout -f ${timeformat}"

export timexec=''

declare -r cfg_plinkctime=$( cfgvar_get plinkctime )
declare -r cfg_plinkmem=$( cfgvar_get plinkmem )

run() {
  local cmd="$*"
  local plinkflag="--memory 100"
  local sbatchflag="--mem-per-cpu=120MB"
  if [ "${cmd}" != "" ] ; then
    declare -r cfg_execommand=$( cfgvar_get execommand )
    case ${cfg_execommand} in
      *bash*)
        echo "running '$( basename ${cmd} )'.."
        export plinkexec="${BASEDIR}/lib/3rd/plink"
        ${timexec} ${cfg_execommand} ${cmd}
        ;;
      *sbatch*)
        echo "submitting '$( basename ${cmd} )'.."
        if [ "${cfg_plinkmem}" != "" ] ; then
          plinkflag="--memory ${cfg_plinkmem}"
          sbatchflag="--mem-per-cpu=$((cfg_plinkmem + cfg_plinkmem/10))MB"
        fi
        export plinkexec="${BASEDIR}/lib/3rd/plink ${plinkflag}"
        jobnum=$( ${cfg_execommand} ${sbatchflag} --time=${cfg_plinkctime} ${jobdepflag} \
                  --output=${opt_outprefix}_slurm_%x_%j.out --parsable ${cmd} )
        joblist=${joblist}:${jobnum}
        jobdepflag="--dependency=afterok${joblist}"
        echo submitted job ${jobnum}.
        ;;
      *)
        echo "error: unknown execution command '${cfg_execommand}'" >&2
        exit 1
        ;;
    esac
  fi
}

#---------------------------------------------------------------------------------------------------

# initialize log file

> "${opt_outprefixbase}.log"

{ 
  echo
  echo -e "================================================================================"
  echo -e "$( basename $0 ) -- $(date)"
  echo -e "================================================================================"
  echo -e "\ngenotype files:\n\n$( printf '  %s\n' ${opt_inputfiles} )\n"
  echo -e "================================================================================"
  echo
} | printlog 1

# print variables
cfgvar_show_config | printlog 1

echo -e "\n================================================================================\n" \
  | printlog 1

#---------------------------------------------------------------------------------------------------

# source utility functions

source "${BASEDIR}/progs/qc-tools.sh"

#---------------------------------------------------------------------------------------------------
# pre-processing
#---------------------------------------------------------------------------------------------------

# export vars

declare  cfg_uid
         cfg_uid="$( cfgvar_get uid )"
readonly cfg_uid
export cfg_uid

declare  cfg_refprefix
         cfg_refprefix="$( cfgvar_get refprefix )"
readonly cfg_refprefix
export cfg_refprefix

export opt_dryimpute
export opt_varwhitelist="${opt_outprefixbase}/.i/qc/a_vwlist.mrk"
export opt_refprefix="${opt_outprefixbase}/.i/qc/a_refset"

#---------------------------------------------------------------------------------------------------

export opt_outprefix="${opt_outprefixbase}/.i/qc/a_vwlist"

# call wlist
${opt_imputeonly} || run "${BASEDIR}/progs/qc-wlist.sh"

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

export opt_outprefix="${opt_outprefixbase}/.i/qc/a_recode"

# call recode
${opt_imputeonly} || run "${BASEDIR}/progs/qc-recode.sh"

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

export opt_inprefix="${opt_outprefixbase}/.i/qc/a_recode"
export opt_outprefix="${opt_outprefixbase}/.i/qc/b_align"

# call align
${opt_imputeonly} || run "${BASEDIR}/progs/qc-align.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# merge batches (if more than a single one present)

export opt_inprefix="${opt_outprefixbase}/.i/qc/b_align"
export opt_outprefix="${opt_outprefixbase}/.i/qc/c_proc"

# call merge
${opt_imputeonly} || run "${BASEDIR}/progs/qc-merge.sh"

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
${opt_imputeonly} || run "${BASEDIR}/progs/init-bio.sh"

# cleanup
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/c_proc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/d_hqset"

# call hqset
${opt_imputeonly} || run "${BASEDIR}/progs/qc-hqset.sh"

# cleanup
unset opt_outprefix
unset opt_inprefix

#---------------------------------------------------------------------------------------------------

# identify duplicate, mixup and related individuals

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/c_proc"
export opt_hqprefix="${opt_outprefixbase}/.i/qc/d_hqset"
export opt_outprefix="${opt_outprefixbase}/.i/qc/e_indqc"
export opt_biofile="${opt_outprefixbase}/bio.txt"

# call indqc
${opt_imputeonly} || run "${BASEDIR}/progs/qc-ind.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

# export vars
export opt_hqprefix="${opt_outprefixbase}/.i/qc/d_hqset"
export opt_inprefix="${opt_outprefixbase}/.i/qc/e_indqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/f_varqc"

# call qcvar
${opt_imputeonly} || run "${BASEDIR}/progs/qc-var.sh"

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
${opt_imputeonly} || run "${BASEDIR}/progs/qc-final.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_batchoutprefix

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/h_hqset"

# call hqset
${opt_imputeonly} || run "${BASEDIR}/progs/qc-hqset.sh" "${cfg_refprefix}"

# cleanup
unset opt_outprefix
unset opt_inprefix

#---------------------------------------------------------------------------------------------------

# compute genetic PCs

# export vars
export opt_hqprefix="${opt_outprefixbase}/.i/qc/h_hqset"
export opt_inprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_biofile="${opt_outprefixbase}/bio.txt"

# call getpcs
${opt_imputeonly} || run "${BASEDIR}/progs/qc-getpcs.sh"

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

${opt_qconly} && { echo "done."; exit 0; }

#---------------------------------------------------------------------------------------------------

# convert

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/g_finqc"
export opt_outprefix="${opt_outprefixbase}/.i/qc/i_pre"

# call convert
${opt_imputeonly} || run "${BASEDIR}/progs/convert.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# impute

# export vars
export opt_inprefix="${opt_outprefixbase}/.i/qc/i_pre"
export opt_outprefix="${opt_outprefixbase}/.i/im/i_imp"

# call impute
run "${BASEDIR}/progs/impute.sh"

# cleanup
unset opt_inprefix
unset opt_outprefix

echo "done."

