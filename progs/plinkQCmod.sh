#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

#-------------------------------------------------------------------------------

# define default configuration

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source ${BASEDIR}/progs/cfgmgr.sh

cfgvar_init_from_file ${BASEDIR}/progs/config.default

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
  [see default configuration file '${BASEDIR}/progs/config.default' for more information]

EOF
exit 0
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
  *)
    usage
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
fi

#---------------------------------------------------------------------------------------------------

# set user configuration, if any, and print all information

[ -s "${opt_cfgfile}" ] && cfgvar_update_from_file ${opt_cfgfile}

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

awk --version 2> /dev/null || { echo 'awk is not installed. aborting..'; exit 0; }
echo
join --version 2> /dev/null || { echo 'join is not installed. aborting..'; exit 0; }
echo
bedtools --version 2> /dev/null || { echo 'bedtools is not installed. aborting..'; exit 0; }
echo
plink --version 2> /dev/null || {
  echo 'plink is not installed.';
  echo "plink is required by $( basename $0 ). plink source codes and builds can be found at"
  echo "www.cog-genomics.org. note that some of the functionalities needed by $( basename $0 )"
  echo "were not implemented in plink2 at the time of writing."
  exit 0
}
echo
R --version 2> /dev/null || { echo 'R is not installed. aborting..'; exit 0; }

echo -e "================================================================================\n"

#---------------------------------------------------------------------------------------------------

# alignment

# export vars
export opt_inputfiles
export opt_minivarset
export opt_samplewhitelist
export opt_batchoutprefix=${opt_outprefixbase}_a_processed_batch
export opt_outprefix=${opt_outprefixbase}_a_align

# call align
bash ${BASEDIR}/progs/align.sh

#---------------------------------------------------------------------------------------------------

# biography

# initialize sample biography file
declare -r cfg_uid="$( cfgvar_get uid )"
cut -f 1,2 ${opt_outprefix}.fam | awk -v uid=${cfg_uid} '
BEGIN{ OFS="\t"; print( uid, "FID", "IID" ) } { print( $1"_"$2, $0 ) }
' | sort -u -k 1,1 > ${opt_outprefixbase}.bio

# ...

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix=${opt_outprefixbase}_a_align
export opt_outprefix=${opt_outprefixbase}_a_hqset

# call gethqset
bash ${BASEDIR}/progs/gethqset.sh

#---------------------------------------------------------------------------------------------------

# identify duplicate and mixup individuals

# export vars
export opt_inprefix=${opt_outprefixbase}_a_align
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_b_clean
export opt_biofile=${opt_outprefixbase}.bio

# call iddupmix.sh
bash ${BASEDIR}/progs/iddupmix.sh

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

# export vars
export opt_inprefix=${opt_outprefixbase}_b_clean
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_c_varqc

# call qcvar.sh
bash ${BASEDIR}/progs/qcvar.sh

#---------------------------------------------------------------------------------------------------

# perform standard individual-level QC

# export vars
export opt_inprefix=${opt_outprefixbase}_c_varqc
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_d_indqc
export opt_biofile=${opt_outprefixbase}.bio

# call qcind.sh
bash ${BASEDIR}/progs/qcind.sh

#---------------------------------------------------------------------------------------------------

# perform final QC (control-HWE? batch effects?)

# export vars
export opt_inprefix=${opt_outprefixbase}_d_indqc
export opt_hqprefix=${opt_outprefixbase}_a_hqset
export opt_outprefix=${opt_outprefixbase}_e_finqc
export opt_batchoutprefix=${opt_outprefixbase}_a_processed_batch

# call qcfinal.sh
bash ${BASEDIR}/progs/qcfinal.sh

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix=${opt_outprefixbase}_e_finqc
export opt_outprefix=${opt_outprefixbase}_e_hqset

# call gethqset.sh
bash ${BASEDIR}/progs/gethqset.sh

#---------------------------------------------------------------------------------------------------

# compute genetic PCs

# export vars
export opt_hqprefix=${opt_outprefixbase}_e_hqset
export opt_biofile=${opt_outprefixbase}.bio

# call getpcs.sh
bash ${BASEDIR}/progs/getpcs.sh

#---------------------------------------------------------------------------------------------------

echo -e "\nall done. check your output files out."
echo -e "\n================================================================================\n"

