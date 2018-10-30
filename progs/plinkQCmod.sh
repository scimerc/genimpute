#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

#---------------------------------------------------------------------------------------------------

# define defaults

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

source ${BASEDIR}/progs/cfgmgr.sh

cfgvar_init_from_file ${BASEDIR}/progs/config.default

#---------------------------------------------------------------------------------------------------

# set user options


#---------------------------------------------------------------------------------------------------

# define vars

declare -r opt_inputfiles="$(cat << EOF
/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/asdtt.bed
/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/qwe.bed
EOF
)"

declare -r opt_mini=1
declare -r opt_refallelesfn="/cluster/projects/p33/nobackup/tmp/refalleles.txt"
declare -r opt_samplewhitelist=""

#---------------------------------------------------------------------------------------------------

# alignment

# export vars
export opt_mini
export opt_refallelesfn
export opt_samplewhitelist
export opt_inputfiles
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_a_align
export opt_batchoutprefix=/cluster/projects/p33/nobackup/tmp/test_a_filtered_batches

# call align
bash ${BASEDIR}/progs/align.sh

#---------------------------------------------------------------------------------------------------

# biography

# initialize sample biography file
declare -r opt_biofile=/cluster/projects/p33/nobackup/tmp/test_0.bio
declare -r uid="000UID"
export uid
cut -f 1,2 ${opt_outprefix}.fam | awk -v uid=${uid} 'BEGIN{
  OFS="\t"; print( uid, "FID", "IID" )
} { print( $1"_"$2, $0 ) }' | sort -u -k 1,1 > ${opt_biofile}

# ...

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix=/cluster/projects/p33/nobackup/tmp/test_a_align
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_a_hqset

# call gethqset
bash ${BASEDIR}/progs/gethqset.sh

#---------------------------------------------------------------------------------------------------

# identify duplicate and mixup individuals

# export vars
export opt_inprefix=/cluster/projects/p33/nobackup/tmp/test_a_align
export opt_hqprefix=/cluster/projects/p33/nobackup/tmp/test_a_hqset
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_b_clean
export opt_biofile

# call iddupmix.sh
bash ${BASEDIR}/progs/iddupmix.sh

#---------------------------------------------------------------------------------------------------

# perform standard variant-level QC

# export vars
export opt_inprefix=/cluster/projects/p33/nobackup/tmp/test_b_clean
export opt_hqprefix=/cluster/projects/p33/nobackup/tmp/test_a_hqset
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_c_varqc

# call qcvar.sh
bash ${BASEDIR}/progs/qcvar.sh

#---------------------------------------------------------------------------------------------------

# perform standard individual-level QC

# export vars
export opt_inprefix=/cluster/projects/p33/nobackup/tmp/test_c_varqc
export opt_hqprefix=/cluster/projects/p33/nobackup/tmp/test_a_hqset
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_d_indqc
export opt_biofile

# call qcind.sh
bash ${BASEDIR}/progs/qcind.sh

#---------------------------------------------------------------------------------------------------

# perform final QC (control HWE? batch effects?)

# export vars
export opt_inprefix=/cluster/projects/p33/nobackup/tmp/test_d_indqc
export opt_hqprefix=/cluster/projects/p33/nobackup/tmp/test_a_hqset
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_e_finqc
export opt_batchoutprefix=/cluster/projects/p33/nobackup/tmp/test_a_filtered_batches

# call qcfinal.sh
bash ${BASEDIR}/progs/qcfinal.sh

#---------------------------------------------------------------------------------------------------

# get high quality set

# export vars
export opt_inprefix=/cluster/projects/p33/nobackup/tmp/test_e_finqc
export opt_outprefix=/cluster/projects/p33/nobackup/tmp/test_e_hqset

# call gethqset.sh
bash ${BASEDIR}/progs/gethqset.sh

#---------------------------------------------------------------------------------------------------

# compute genetic PCs

# export vars
export opt_hqprefix=/cluster/projects/p33/nobackup/tmp/test_e_finqc

# call getpcs.sh
bash ${BASEDIR}/progs/getpcs.sh

#---------------------------------------------------------------------------------------------------


exit 0


# make directory to accommodate output prefix
# mydir=$( dirname ${myprefix} )
# mkdir -p ${mydir}


AWKPATH="${AWKPATH}:${BASEDIR}/lib/awk"
AWKLOCINCLUDE=$( printf -- '-f %s\n' $( echo ${AWKINCLUDE} \
  abs.awk \
  nucleocode.awk \
  genotype.awk \
  gflip.awk \
  gmatch.awk \
  round.awk \
) | sort -u )
alias awk='awk --lint=true'
alias join='join --check-order'
shopt -s expand_aliases
fdrscript="R --slave --file=${BASEDIR}/progs/fdr.Rscript"
hvmscript="R --slave --file=${BASEDIR}/progs/het_vs_miss.Rscript"
plinkexe=$( ( which plink || true ) 2> /dev/null )
if [[ "${plinkexe}" == "" ]] ; then
  echo "plink is required by $( basename $0 ). plink source codes and builds can be found at"
  echo "www.cog-genomics.org. note that some of the functionalities needed by $( basename $0 )"
  echo "were not implemented in plink2 at the time of writing."
  exit
fi

opt_parser="${BASEDIR}/lib/sh/opt_parser.bash"
opt_list=("b=" "dups" "g=" "i=" "maf=" "nohvm" "o=" "p=" "r=" "s=" "v=" "x" "h")
get_opt ()
{
  case $1 in
    "b" )
      blacklist=`echo "$2" | sed "s/^=\+//"` ;;
    "dups" )
      keepdups='on' ;;
    "g" )
      genome=`echo "$2" | sed "s/^=\+//"` ;;
    "i" )
      samplefile=`echo "$2" | sed "s/^=\+//"` ;;
    "maf" )
      myfreq_std=`echo "$2" | sed "s/^=\+//"` ;;
    "nohvm" )
      hvm='off' ;;
    "o" )
      myprefix=`echo "$2" | sed "s/^=\+//"` ;;
    "p" )
      phenotypes=`echo "$2" | sed "s/^=\+//"` ;;
    "r" )
      refalleles=`echo "$2" | sed "s/^=\+//"` ;;
    "s" )
      samplemiss=`echo "$2" | sed "s/^=\+//"` ;;
    "v" )
      varmiss=`echo "$2" | sed "s/^=\+//"` ;;
    "x" )
      mini="yes" ;;
    "h" )
      helpme='yes' ;;
  esac
}
helpme=''
blacklist=''
genome='b37'
hvm='on'
hweneglogp=12
hweneglogp_ctrl=4
hwflag='midp include-nonctrl'
keepdups='off'
minindcount=100
minvarcount=100
mini='no'
myfreq_hq=0.2
myfreq_std=0.05
myprefix='plink'
phenotypes=''
pihat=0.9
refalleles=''
samplefile=''
samplemiss=0.05
varmiss=0.05
tmpargs=''
while [[ "${tmpargs}" == "" ]] ; do
  tmpargs=$( mktemp .tmpXXXXXXXX )
done
source ${opt_parser} > ${tmpargs} || true
mybatches=$( cat ${tmpargs} )
rm -f ${tmpargs}

if [[ "${mybatches}" == "" || "${helpme}" != "" ]] ; then

  if [[ "${helpme}" == "" ]] ; then echo "missing input."; fi

  echo -e "\n USAGE:"
  echo -e "   $( basename $0 ) [OPTIONS] <bed|bcf|vcf file(s)>\n"
  echo -e "   where\n"
  echo -e "   <bed|bcf|vcf file(s)> are the genotype files to be merged."
  echo -e "\n OPTIONS:"
  echo -e "   -b <black list>      recommended list of bad (high LD) regions (bed format)"
  echo -e "   -dups            keep duplicate individuals [default: off]"
  echo -e "   -g <genome version>    genome version [default: hg19]"
  echo -e "   -i <sample file>       optional individual selection file"
  echo -e "   -maf <freq>        minor allele frequency filter"
  echo -e "   -nohvm           turns off heterozygosity VS missingness tests"
  echo -e "   -o <output prefix>     optional output prefix [default: 'plink']"
  echo -e "   -p <phenotype file>    optional file with 'affected'/'control' status tags"
  echo -e "   -r <reference alleles>   tab separated reference alleles"
  echo -e "                [format: <CHR:POS> <A1> <A2>]"
  echo -e "   -s <sample missingness>  max. fraction of variants missing in an individual"
  echo -e "   -v <variant missingness>   max. fraction of genotypes missing for one variant"
  echo -e "   -h             print help\n"

  exit

fi

module load R

echo -e "==== preMaCH.sh -- $(date)\n"
echo -e "==== options in effect:\n"
echo "   blacklist=${blacklist}"
echo "   keepdups=${keepdups}"
echo "   genome=${genome}"
echo "   samplefile=${samplefile}"
echo "   hvm=${hvm}"
echo "   myprefix=${myprefix}"
echo "   phenotypes=${phenotypes}"
echo "   refalleles=${refalleles}"
echo "   samplemiss=${samplemiss}"
echo "   varmiss=${varmiss}"
echo -e "\n==================================\n"
echo -e "batch files:\n$( ls ${mybatches} )"
echo -e "\n==================================\n"

awk --version
echo
join --version
echo
R --version






echo -e "\nall done. check your output files out."
echo -e "\n================================================================================\n"

