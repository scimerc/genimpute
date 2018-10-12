#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

#-------------------------------------------------------------------------------

# define defaults

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
export BASEDIR

#-------------------------------------------------------------------------------

# set user options


#-------------------------------------------------------------------------------

# define vars

declare -r opt_inputfiles="$(cat << EOF
/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/asdt.bed
/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/qwe.bed
EOF
)"

# set hardy weinberg test p-value threshold for the general case
# if [[ "${phenotypes}" == "" || ! -f "${phenotypes}" ]] ; then
#   hweneglogp=${hweneglogp_ctrl}
# fi

declare -r opt_mini=1
declare -r opt_refallelesfn=""
declare -r opt_samplewhitelist=""

#-------------------------------------------------------------------------------

# alignment

# export vars

export opt_mini
export opt_refallelesfn
export opt_samplewhitelist
export opt_inputfiles
export opt_outprefix=/tmp/test_aligned
export opt_batchoutprefix=/tmp/test_filtered_batches

# call align
bash ${BASEDIR}/progs/align.sh

#-------------------------------------------------------------------------------

# biography

# initialize sample biography file
declare -r opt_biofile=/tmp/test.bio
uid="000UID"
cut -f 1,2 ${opt_outprefix}.fam | awk -v uid=${uid} 'BEGIN{
  OFS="\t"; print( uid, "FID", "IID" )
} { print( $1"_"$2, $0 ) }' | sort -u -k 1,1 > ${opt_biofile}

# ...

#-------------------------------------------------------------------------------
# get high quality set

# export vars

export opt_inprefix=/tmp/test_aligned
export opt_outprefix=/tmp/test_hqset
export opt_biofile

# call gethqset
bash ${BASEDIR}/progs/gethqset.sh

#-------------------------------------------------------------------------------


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

# define filenames
# input: merged plink set
# output: hq plink set (with imputed sex, if possible)
# 1) get sex hq-variants from input file
# 2) get non-sex hq-variants from input file
# 3) extract all hq variants from input file and make hq plink set
# 4) LD-prune hq variants from 3)
# 5) impute sex once with all standard hq variants from 4
# 6) if sex could be imputed for enough individuals, then
#      impute it once again after HWE tests
# 7) update biography file with sex information



# XXX) if sex could be imputed for enough individuals, then
#      set --remove plink flag for sex-undetermined individuals
#      set --update-sex plink flag
#     usex="--update-sex ${prunehqprefix}_isex.fam 3"
#     awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' \
#       ${prunehqprefix}_isex.fam > ${prunehqprefix}_isex.nosex
#     if [[ -s "${prunehqprefix}_isex.nosex" ]] ; then
#       nosex="--remove ${prunehqprefix}_isex.nosex"
#     fi


keepflag=''
echo "checking sample quality.."
if [[ ! -s "${uniquefile}" ]] ; then
  cp ${prunehqprefix}.fam ${cleanfile}
  if [[ "${hvm}" == "on" ]] ; then
    echo "computing individual heterozygosity and missing rates.."
    plink --bfile ${prunehqprefix} --het --out ${hqprefix}_sq
    plink --bfile ${prunehqprefix} --missing --out ${hqprefix}_sq
    ${hvmscript} --args -m ${hqprefix}_sq.imiss -h ${hqprefix}_sq.het -o ${hqprefix}
    tmpbiofile=$( mktemp ${mybiofile}.tmpXXXX )
    (
      cut -f 3 ${cleanfile} | sort -u | join -t $'\t' -v1 ${mybiofile} - | awk -F $'\t' '{
        OFS="\t"
        if ( NR == 1 ) print( $0, "het_VS_miss" )
        else print( $0, 0 )
      }'
      cut -f 3 ${cleanfile} | sort -u | join -t $'\t' ${mybiofile} - | awk -F $'\t' '{
        OFS="\t"
        print( $0, 1 )
      }'
    ) | sort -t $'\t' -u -k 1,1 > ${tmpbiofile}
    mv ${tmpbiofile} ${mybiofile}
  fi
  if [[ -s "${cleanfile}" ]] ; then
    genomefile=${myprefix}.genome.gz
    echo "identifying duplicate individuals.."
    if [[ "${keepdups}" == "on" ]] ; then pihat=1 ; fi
    plink --bfile ${prunehqprefix} --keep ${cleanfile} --genome gz --out ${myprefix}
    plink --bfile ${prunehqprefix} --keep ${cleanfile} --cluster --read-genome ${genomefile} \
        --rel-cutoff ${pihat} --out ${hqprefixunique}
    tmpbiofile=$( mktemp ${mybiofile}.tmpXXXX )
    zcat ${genomefile} | sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;' \
    | awk -F $'\t' -v uid=${uid} 'BEGIN{
      OFS="\t"
      printf( "%s\tRSHIP\n", uid )
    } {
      if ( NR>1 && $10>=0.1 ) {
        uid0=$1"_"$2
        uid1=$3"_"$4
        if ( uid0 in relarr )
          relarr[uid0] = relarr[uid0]","uid1"("$10")"
        else relarr[uid0] = uid1"("$10")"
        if ( uid1 in relarr )
          relarr[uid1] = relarr[uid1]","uid0"("$10")"
        else relarr[uid1] = uid0"("$10")"
      }
    } END{
      for ( uid in relarr )
        print( uid, relarr[uid] )
    }' | sort -t $'\t' -u -k 1,1 | join -t $'\t' -a1 -e '-' ${mybiofile} - > ${tmpbiofile}
    mv ${tmpbiofile} ${mybiofile}
  fi
fi
if [[ -s "${uniquefile}" ]] ; then
  keepflag="--keep ${uniquefile}"
  tmpbiofile=$( mktemp ${mybiofile}.tmpXXXX )
  (
    awk -F $'\t' '{ print( $1"\t"$2 ); }' ${uniquefile} | sort -u \
    | join -t $'\t' -v1 ${mybiofile} - | awk -F $'\t' '{
      OFS="\t"
      if ( NR == 1 ) print( $0, "duplicate" )
      else print( $0, 1 )
    }'
    awk -F $'\t' '{ print( $1"\t"$2 ); }' ${uniquefile} | sort -u \
    | join -t $'\t' ${mybiofile} - | awk -F $'\t' '{
      OFS="\t"
      print( $0, 1 )
    }'
  ) | sort -t $'\t' -u -k 1,1 > ${tmpbiofile}
  mv ${tmpbiofile} ${mybiofile}
fi
outprefix=${myprefix}_clean
if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
  plink --bfile ${myprefix} ${usex} ${keepflag} --make-bed --out ${outprefix}
  if [[ -f "${outprefix}.bim" ]] ; then
    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
  fi
  if [[ -f "${outprefix}.fam" ]] ; then
    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
  fi
fi
myprefix=${outprefix}

outprefix=${myprefix}_varQC
if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
  tmpvarmiss=${varmiss}
  N=$( wc -l ${myprefix}.fam | cut -d ' ' -f 1 )
  if (( N < minindcount )) ; then tmpvarmiss=0.1 ; fi
  hweneglogp_sex=$(( hweneglogp*2 ))
  if (( hweneglogp_sex > 12 )) ; then hweneglogp_sex=12 ; fi
  plink --bfile ${myprefix} ${nosex} --not-chr 23,24 --geno ${tmpvarmiss} --maf ${myfreq_std} \
    --hwe 1.E-${hweneglogp} ${hwflag} --make-just-bim --out ${outprefix}_nonsex
  plink --bfile ${myprefix} ${nosex} --chr 23,24 --geno ${tmpvarmiss} --maf ${myfreq_std} \
    --hwe 1.E-${hweneglogp_sex} ${hwflag} --make-just-bim --out ${outprefix}_sex
  if [[ -s "${outprefix}_nonsex.bim" || -s "${outprefix}_sex.bim" ]] ; then
    cut -f 2 ${outprefix}_*sex.bim | sort -u > ${outprefix}.mrk
  fi
  if [[ -s "${outprefix}.mrk" ]] ; then
    plink --bfile ${myprefix} --extract ${outprefix}.mrk --make-bed --out ${outprefix}
    if [[ -f "${outprefix}.bim" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
    fi
    if [[ -f "${outprefix}.fam" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
    fi
  fi
fi
myprefix=${outprefix}

outprefix=${myprefix}_sampleQC
if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
  tmpsamplemiss=${samplemiss}
  N=$( wc -l ${myprefix}.bim | cut -d ' ' -f 1 )
  if (( N < minindcount )) ; then tmpsamplemiss=0.1 ; fi
  plink --bfile ${myprefix} --mind ${tmpsamplemiss} --make-bed --out ${outprefix}
  if [[ -f "${outprefix}.bim" ]] ; then
    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
  fi
  if [[ -f "${outprefix}.fam" ]] ; then
    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
  fi
  tmpbiofile=$( mktemp ${mybiofile}.tmpXXXX )
  (
    awk -F $'\t' '{ print( $1"\t"$2 ); }' ${outprefix}.fam | sort -u \
    | join -t $'\t' -v1 ${mybiofile} - | awk -F $'\t' '{
      OFS="\t"
      if ( NR == 1 ) print( $0, "coverage" )
      else print( $0, 0 )
    }'
    awk -F $'\t' '{ print( $1"\t"$2 ); }' ${outprefix}.fam | sort -u \
    | join -t $'\t' ${mybiofile} - | awk -F $'\t' '{
      OFS="\t"
      print( $0, 1 )
    }'
  ) | sort -t $'\t' -u -k 1,1 > ${tmpbiofile}
  mv ${tmpbiofile} ${mybiofile}
fi
myprefix=${outprefix}

if [[ "${phenotypes}" != "" && -f "${phenotypes}" ]] ; then
  outprefix=${myprefix}_pheno
  if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
    plink --bfile ${myprefix} --make-pheno ${phenotypes} affected --make-bed \
      --out ${outprefix}
    if [[ -f "${outprefix}.bim" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
    fi
    if [[ -f "${outprefix}.fam" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
    fi
  fi
  myprefix=${outprefix}
  outprefix=${myprefix}_hwe_ctrl
  if [[ ! -f "${outprefix}.bed" || ! -f "${outprefix}.bim" || ! -f "${outprefix}.fam" ]] ; then
    hweneglogp_ctrl_sex=$(( hweneglogp_ctrl/2 ))
    if (( hweneglogp_ctrl_sex < 12 )) ; then hweneglogp_ctrl_sex=12 ; fi
    plink --bfile ${myprefix} ${nosex} --not-chr 23,24 --hwe 1.E-${hweneglogp_ctrl} \
      --make-just-bim --out ${outprefix}_nonsex
    plink --bfile ${myprefix} ${nosex} --chr 23,24 --hwe 1.E-${hweneglogp_ctrl_sex} \
      --make-just-bim --out ${outprefix}_sex
    if [[ -s "${outprefix}_nonsex.bim" || -s "${outprefix}_sex.bim" ]] ; then
      cut -f 2 ${outprefix}_*sex.bim | sort -u > ${outprefix}.mrk
    fi
    if [[ -s "${outprefix}.mrk" ]] ; then
      plink --bfile ${myprefix} --extract ${outprefix}.mrk --make-bed --out ${outprefix}
      if [[ -f "${outprefix}.bim" ]] ; then
        perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
      fi
      if [[ -f "${outprefix}.fam" ]] ; then
        perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
      fi
    fi
  fi
  myprefix=${outprefix}
  mycontrols=${mydir}/$( basename ${phenotypes} .txt )_ctrl.txt
  mycontrols_with_batch=${mydir}/$( basename ${phenotypes} .txt )_ctrl.batch
  awk -F $'\t' '{ if ( tolower($3) == "control" ) print; }' ${phenotypes} > ${mycontrols}
else
  mycontrols=${myprefix}_ctrl.txt
  mycontrols_with_batch=${myprefix}_ctrl.batch
  cut -f 1,2,6 ${myprefix}.fam | perl -p -e 's/[ \t]+/\t/g' > ${mycontrols}
fi

if (( ${#mybatchvec[*]} > 1 )) ; then
  tmpfile=$( mktemp .tmpbcXXXXX )
  echo "created temporary batch control file ${tmpfile}."
  echo -n > ${mycontrols_with_batch}
  for batch in $( cat ${mybatchfile} ) ; do
    tmpbatch=${mydir}/$( basename ${batch} )_ctrl
    plink --bfile ${batch} --keep ${mycontrols} --make-bed --out ${tmpbatch}
    if [[ -f "${tmpbatch}.bim" && -f "${tmpbatch}.fam" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${tmpbatch}.bim ${tmpbatch}.fam
      awk -F $'\t' -v batch=$( basename ${batch} ) '{ OFS="\t"; print( $1, $2, batch ) }' \
        ${tmpbatch}.fam >> ${mycontrols_with_batch}
    fi
  done
  sort -u -k 1,2 ${mycontrols_with_batch} > ${tmpfile}
  mv ${tmpfile} ${mycontrols_with_batch}
  tmpfile=$( mktemp .tmpmodXXXXX )
  echo "created temporary model file ${tmpfile}."
  bevarfile=${myprefix}.exclude
  echo -n > ${bevarfile}
  for bbatch in $( cat ${mybatchfile} ) ; do
    batch=$( basename ${bbatch} )
    echo "assessing batch effects for '${batch}'.."
    plinkfile=${mydir}/plink_${batch}
    plink --bfile ${myprefix} --make-pheno ${mycontrols_with_batch} ${batch} --model \
      --out ${plinkfile}
    if [[ -f "${plinkfile}.model" ]] ; then
      perl -p -e 's/^[ \t]+//g;' ${plinkfile}.model | perl -p -e 's/[ \t]+/\t/g' > ${tmpfile}
      mv ${tmpfile} ${plinkfile}.model
      for atest in $( cut -f 5 ${plinkfile}.model | tail -n +2 | sort -u ) ; do
        echo "summarizing ${atest} tests.."
        awk -F $'\t' -v atest=${atest} 'NR > 1 && $5 == atest' ${plinkfile}.model \
        | cut -f 2,10 | ${fdrscript} > ${tmpfile}
        awk -F $'\t' '{ if ( $2 < 0.9 ) print( $1 ); }' ${tmpfile} >> ${bevarfile}
      done
    fi
  done
  sort -u ${bevarfile} > ${tmpfile}
  mv ${tmpfile} ${bevarfile}
  outprefix=${myprefix}_nbe
  plink --bfile ${myprefix} --exclude ${bevarfile} --make-bed --out ${outprefix}
  if [[ -f "${outprefix}.bim" ]] ; then
    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
  fi
  if [[ -f "${outprefix}.fam" ]] ; then
    perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
  fi
fi
if [[ -f "${outprefix}.fam" ]] ; then
  cp ${outprefix}.fam ${outprefix}.fam.org
  awk '{ OFS="\t"; for ( k = 1; k < 5; k++ ) gsub( "[/+-]", "_", $k ); print; }' \
  ${outprefix}.fam.org > ${outprefix}.fam
fi

hqprefix=${outprefix}_hq
genomefile=${hqprefix}.genome.gz
plink --bfile ${outprefix} ${excludeopt} --maf ${myfreq_hq} --make-bed --out ${hqprefix}
plink --bfile ${hqprefix} --indep-pairphase 500 5 0.2 --out ${hqprefix}_LD
plink --bfile ${hqprefix} --extract ${hqprefix}_LD.prune.in --make-bed --out ${hqprefix}_LDpruned
plink --bfile ${hqprefix}_LDpruned --genome gz --out ${hqprefix}
plink --bfile ${hqprefix}_LDpruned --cluster --read-genome ${genomefile} \
    --pca header tabs var-wts --out ${hqprefix}

echo -e "\nall done. check your output files out."
echo -e "\n================================================================================\n"

