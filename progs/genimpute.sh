#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

#-------------------------------------------------------------------------------

# define default configuration

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source ${BASEDIR}/progs/cfgmgr.sh

cfgvar_init_from_file ${BASEDIR}/progs/cfgqc.default

#---------------------------------------------------------------------------------------------------

# define default options and parse command line

declare opt_dryimpute=0
declare opt_minivarset=0
declare opt_outprefixdefault='plinkqc'
declare opt_outprefixbase=${opt_outprefixdefault}
declare opt_samplewhitelist=""

usage() {
cat << EOF

USAGE: $( basename $0 ) [OPTIONS] <bed|bcf|vcf file(s)>

  where <bed|bcf|vcf file(s)> are the genotype files to be merged and qc'd.


OPTIONS:
  -c <config file>      optional configuration file
  -d                    dry imputation: write scripts but do not run them
  -m                    reduce variant set to the minimal one common to all
  -w <sample file>      optional white list of individuals to restrict qc to
  -o <output prefix>    optional output prefix [default: '${opt_outprefixdefault}']
  -h                    show help message


CONFIGURATION:
  the <config file> may contain custom definitions for a number of variables.
  [see default configuration file 'cfgqc.default' for more information]

EOF
}

debugout() {
  local -r lvl
  lvl=$1
  local -r lvl_max
  lvl_max ="$( cfgvar_get debug_lvl )"
  # print struff to log
  while read line; do
    # check log-level
    if [ "$lvl" -le "$lvl_max" ]; then
      echo $(date) "$line"
    fi
  done
  return 0
}
export -f debugout

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

# export input and output names
export opt_inputfiles
export opt_outprefixbase

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

echo
echo -e "================================================================================"
echo -e "$( basename $0 ) -- $(date)"
echo -e "================================================================================"
echo -e "\ngenotype files:\n$( ls -1 ${opt_inputfiles} )\n"
echo -e "================================================================================"
echo

# lock variables and print them
cfgvar_setall_readonly
cfgvar_show_config
echo

echo -e "================================================================================\n"

#---------------------------------------------------------------------------------------------------

# define executables

plinkexec="${BASEDIR}/lib/3rd/plink --allow-extra-chr"
export plinkexec

#---------------------------------------------------------------------------------------------------

# check dependencies

locale
echo

awk --version 2> /dev/null || { echo 'awk is not installed. aborting..'; exit 1; }
echo
join --version 2> /dev/null || { echo 'join is not installed. aborting..'; exit 1; }
echo
R --version 2> /dev/null || { echo 'R is not installed. aborting..'; exit 1; }

echo -e "================================================================================\n"

#---------------------------------------------------------------------------------------------------

# source utility functions

source ${BASEDIR}/progs/qc-tools.sh

#---------------------------------------------------------------------------------------------------

# pre-processing

# export vars
export opt_dryimpute
export opt_minivarset
export opt_samplewhitelist
export opt_varwhitelist=${opt_outprefixbase}_vwlist.mrk
export opt_outprefix=${opt_outprefixbase}_vwlist
export opt_refprefix=${opt_outprefixbase}_refset

# call wlist?
if [ $opt_minivarset -eq 1 ] ; then
  # intersect batch variant sets to generate the minimal variant whitelist
  bash ${BASEDIR}/progs/qc-wlist.sh
fi

#---------------------------------------------------------------------------------------------------

export opt_outprefix=${opt_outprefixbase}_a_recode

# call recode
bash ${BASEDIR}/progs/qc-recode.sh

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

export opt_inprefix=${opt_outprefixbase}_a_recode
export opt_outprefix=${opt_outprefixbase}_b_align

# call align
bash ${BASEDIR}/progs/qc-align.sh

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

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
  if [ ${Norg} -ne ${Nuid} ] ; then
    echo '========> conflicts generated in universal IDs.'
    echo 'please recode IDs so they do not lead to conflicts.'
    exit 1
  fi
  unset opt_outprefix
  unset opt_biofile
fi

#---------------------------------------------------------------------------------------------------

# get high quality set

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

# export vars
export opt_inprefix=${opt_outprefixbase}_c_proc
export opt_hqprefix=${opt_outprefixbase}_d_hqset
export opt_outprefix=${opt_outprefixbase}_e_indqc
export opt_biofile=${opt_outprefixbase}.bio

# call hqset and mixrel
bash ${BASEDIR}/progs/qc-mixrel.sh

# cleanup
unset opt_hqprefix
unset opt_inprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

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

unset opt_refprefix

echo -e "\nquality control ultimated. you may check your output files out.\n"

#---------------------------------------------------------------------------------------------------

# impute

# export vars
export opt_inprefix=${opt_outprefixbase}_g_finqc
export opt_outprefix=${opt_outprefixbase}_i_imp

# call impute
bash ${BASEDIR}/progs/impute.sh

echo -e "\nrunning imputation..\n"

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

