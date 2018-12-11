#!/usr/bin/env bash

# exit on error
set -Eeou pipefail
trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

#-------------------------------------------------------------------------------

# define default configuration

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source ${BASEDIR}/progs/cfgmgr.sh

cfgvar_init_from_file ${BASEDIR}/progs/cfgqc.default

#---------------------------------------------------------------------------------------------------

# define default options and parse command line

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
  -m                    reduce variant set to the minimal one common to all
  -w <sample file>      optional white list of individuals to restrict qc to
  -o <output prefix>    optional output prefix [default: '${opt_outprefixdefault}']
  -h                    show help message


CONFIGURATION:
  the <config file> may contain custom definitions for a number of variables.
  [see default configuration file 'cfgqc.default' for more information]

EOF
}

while getopts "c:mo:w:h" opt; do
case "${opt}" in
  c)
    opt_cfgfile="${OPTARG}"
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

[ ! -z "${opt_cfgfile+x}" ] && [ -s "${opt_cfgfile}" ] && cfgvar_update_from_file "${opt_cfgfile}"

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

echo -e "\n================================================================================\n"

#---------------------------------------------------------------------------------------------------

# check dependencies

locale
echo

awk --version 2> /dev/null || { echo 'awk is not installed. aborting..'; exit 1; }
echo
join --version 2> /dev/null || { echo 'join is not installed. aborting..'; exit 1; }
echo
plink --version 2> /dev/null || {
  echo 'plink is not installed.';
  echo "plink is required by $( basename $0 ). plink source codes and builds can be found at"
  echo "www.cog-genomics.org. note that some of the functionalities needed by $( basename $0 )"
  echo "were not implemented in plink2 at the time of writing."
  exit 1
}
echo
R --version 2> /dev/null || { echo 'R is not installed. aborting..'; exit 1; }

echo -e "================================================================================\n"

#---------------------------------------------------------------------------------------------------

# source utility functions

source ${BASEDIR}/progs/qc-tools.sh

#---------------------------------------------------------------------------------------------------

# pre-processing

# export vars
export opt_minivarset
export opt_samplewhitelist
export opt_varwhitelist=${opt_outprefixbase}_a_vwlist.mrk
export opt_outprefix=${opt_outprefixbase}_a_intsec

# call wlist?
if [ $opt_minivarset -eq 1 ] ; then
  # intersect batch variant sets to generate the minimal variant whitelist
  bash ${BASEDIR}/progs/qc-wlist.sh
fi

#---------------------------------------------------------------------------------------------------

export opt_outprefix=${opt_outprefixbase}_a_recode

# call recode
bash ${BASEDIR}/progs/qc-recode.sh

#---------------------------------------------------------------------------------------------------

export opt_inprefix=${opt_outprefixbase}_a_recode
export opt_outprefix=${opt_outprefixbase}_a_align

# call align
bash ${BASEDIR}/progs/qc-align.sh

#---------------------------------------------------------------------------------------------------

export opt_inprefix=${opt_outprefixbase}_a_align
export opt_outprefix=${opt_outprefixbase}_a_proc

# call merge
bash ${BASEDIR}/progs/qc-merge.sh

#---------------------------------------------------------------------------------------------------

# cleanup
unset opt_inprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# initialize qc iteration counter
declare qciter=0

# biography

export opt_outprefix=${opt_outprefixbase}_a_proc

# initialize sample biography file
declare -r cfg_uid="$( cfgvar_get uid )"
cut -f 1,2 ${opt_outprefix}.fam | awk -v uid=${cfg_uid} '
BEGIN{ OFS="\t"; print( uid, "FID", "IID" ) } { print( $1"_"$2, $0 ) }
' | sort -u -k 1,1 > ${opt_outprefixbase}.bio

# ...

# cleanup
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# get high quality set and identify duplicate, mixup and related individuals

# export vars
export opt_inprefix=${opt_outprefixbase}_a_proc
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_b_clean
export opt_biofile=${opt_outprefixbase}.bio

# call hqset and mixrel
bash ${BASEDIR}/progs/qc-hqset.sh
bash ${BASEDIR}/progs/qc-mixrel.sh

# cleanup
unset opt_inprefix
unset opt_hqprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

# export vars
export opt_inprefix=${opt_outprefixbase}_b_clean
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_c_varqc

# call qcvar
bash ${BASEDIR}/progs/qc-var.sh

# cleanup
unset opt_inprefix
unset opt_hqprefix
unset opt_outprefix

#---------------------------------------------------------------------------------------------------

# perform standard individual-level QC

# export vars
export opt_inprefix=${opt_outprefixbase}_c_varqc
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_d_indqc
export opt_biofile=${opt_outprefixbase}.bio

# call qcind
bash ${BASEDIR}/progs/qc-ind.sh

# cleanup
unset opt_inprefix
unset opt_hqprefix
unset opt_outprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

# perform final QC (control-HWE? batch effects?)

# export vars
export opt_inprefix=${opt_outprefixbase}_d_indqc
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_e_finqc
export opt_batchoutprefix=${opt_outprefixbase}_a_proc_batch

# call qcfinal
bash ${BASEDIR}/progs/qc-final.sh

# cleanup
unset opt_inprefix
unset opt_hqprefix
unset opt_outprefix
unset opt_batchoutprefix

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix=${opt_outprefixbase}_e_finqc
export opt_hqprefix=${opt_outprefixbase}_e_hqset

# call hqset
bash ${BASEDIR}/progs/qc-hqset.sh

# cleanup
unset opt_inprefix
unset opt_hqprefix

#---------------------------------------------------------------------------------------------------

# compute genetic PCs

# export vars
export opt_hqprefix=${opt_outprefixbase}_e_hqset
export opt_biofile=${opt_outprefixbase}.bio

# call getpcs
bash ${BASEDIR}/progs/qc-getpcs.sh

# cleanup
unset opt_hqprefix
unset opt_biofile

#---------------------------------------------------------------------------------------------------

echo -e "\nall done. check your output files out."
echo -e "\n================================================================================\n"

