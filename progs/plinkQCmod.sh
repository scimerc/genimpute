#!/bin/bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

#-------------------------------------------------------------------------------

# define defaults

# get parent dir of this script
declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"

#-------------------------------------------------------------------------------

# set user options


#-------------------------------------------------------------------------------

# define vars

declare -r inputfiles="$(cat << EOF
/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/asdt.bed
/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/qwe.bed
EOF
)"

declare -r outprefix=/tmp/aligntest

declare -r opt_mini=1
declare -r opt_refallelesfn=""
declare -r opt_samplewhitelist=""

#-------------------------------------------------------------------------------

# export vars

export BASEDIR
export outprefix
export opt_mini
export opt_refallelesfn
export opt_samplewhitelist
export inputfiles

# call align
bash ${BASEDIR}/progs/align.sh

#-------------------------------------------------------------------------------

# export ...

# call ...

#-------------------------------------------------------------------------------

exit 1


# make directory to accommodate output prefix
# mydir=$( dirname ${myprefix} )
# mkdir -p ${mydir}

# master:
# files needed for batch effect qc:
export batchoutprefix=/tmp/...

# align.sh:
plinkoutputfn=${batchoutprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_filtered

# master:
ls ${batchoutprefix}*

# set hardy weinberg test p-value threshold for the general case
# if [[ "${phenotypes}" == "" || ! -f "${phenotypes}" ]] ; then
#   hweneglogp=${hweneglogp_ctrl}
# fi

# initialize sample biography file
# mybiofile=${outprefix}.bio
# uid="000UID"
# cut -f 1,2 ${myprefix}.fam | awk -v uid=${uid} 'BEGIN{
#   OFS="\t"; print( uid, "FID", "IID" )
# } { print( $1"_"$2, $0 ) }' | sort -u -k 1,1 > ${mybiofile}

# pre-process sex chromosomes variants
# parcount=$( awk '$1 == 25' ${myprefix}.bim | wc -l )
# if (( parcount == 0 )) ; then
#   outprefix=${myprefix}_sx
#   plink --bfile ${myprefix} --split-x ${genome} no-fail --make-bed --out ${outprefix}
# fi
# myprefix=${outprefix}


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

mydir=$( dirname ${myprefix} )
mkdir -p ${mydir}

fflag='--bfile'
outprefix=${myprefix}
if (( ${#mybatchvec[*]} > 1 )) ; then
  outprefix=${myprefix}_merged
fi
mybatchfile=${myprefix}_batches
ls ${mybatches} | sed -r 's/[.]bed$//g;' > ${mybatchfile}
mybatchvec=( $( cat ${mybatchfile} ) )

if [[ "${phenotypes}" == "" || ! -f "${phenotypes}" ]] ; then
  hweneglogp=${hweneglogp_ctrl}
fi

myprefix=${outprefix}
mybiofile=${outprefix}.bio
# initialize sample biography file
uid="000UID"
cut -f 1,2 ${myprefix}.fam | awk -v uid=${uid} 'BEGIN{
  OFS="\t"; print( uid, "FID", "IID" )
} { print( $1"_"$2, $0 ) }' | sort -u -k 1,1 > ${mybiofile}
# pre-process sex chromosomes variants
parcount=$( awk '$1 == 25' ${myprefix}.bim | wc -l )
if (( parcount == 0 )) ; then
  outprefix=${myprefix}_sx
  plink --bfile ${myprefix} --split-x ${genome} no-fail --make-bed --out ${outprefix}
fi
myprefix=${outprefix}
hqprefix=${myprefix}_hq
hqprefixunique=${myprefix}_hq_unique
cleanfile=${hqprefix}.clean.id
uniquefile=${hqprefixunique}.rel.id
excludeopt=""
if [[ "${blacklist}" != "" ]] ; then
  excludeopt="--exclude range ${blacklist}"
fi
prunehqprefix=${hqprefix}_LDpruned
if [[ 
    ! -f "${prunehqprefix}.bed" ||
    ! -f "${prunehqprefix}.bim" ||
    ! -f "${prunehqprefix}.fam"
]] ; then
  plink --bfile ${myprefix} --not-chr 23,24 ${excludeopt} --geno ${varmiss} --maf ${myfreq_hq} \
    --hwe 1.E-${hweneglogp_ctrl} ${hwflag} --make-just-bim --out ${hqprefix}_nonsex
  plink --bfile ${myprefix} --chr 23,24 ${excludeopt} --geno ${varmiss} --maf ${myfreq_hq} \
    --make-just-bim --out ${hqprefix}_sex
  if [[ -s "${hqprefix}_nonsex.bim" || -s "${hqprefix}_sex.bim" ]] ; then
    cut -f 2 ${hqprefix}_*sex.bim | sort -u > ${hqprefix}.mrk
    plink --bfile ${myprefix} --extract ${hqprefix}.mrk --make-bed --out ${hqprefix}
  fi
  if [[ -s "${hqprefix}.bed" ]] ; then
    plink --bfile ${hqprefix} --indep-pairphase 500 5 0.2 --out ${hqprefix}_LD
    plink --bfile ${hqprefix} --extract ${hqprefix}_LD.prune.in --make-bed --out ${prunehqprefix}
  fi
fi
usex=''
nosex=''
touch ${prunehqprefix}.bim
xcount=$( awk '$1 == 23' ${prunehqprefix}.bim | wc -l )
if (( xcount > minvarcount )) ; then
  touch ${prunehqprefix}_isex.fam
  # impute sex once with all standard high quality variants
  plink --bfile ${prunehqprefix} --impute-sex --make-bed --out ${prunehqprefix}_isex
  xcount=$( awk '$5 == 1 || $5 == 2' ${prunehqprefix}_isex.fam | wc -l )
  # if sex could be imputed for enough individuals impute it once again after HWE tests
  if (( xcount > minindcount )) ; then
    plink --bfile ${prunehqprefix}_isex --hwe 1.E-${hweneglogp_ctrl} ${hwflag} \
      --make-just-bim --out ${prunehqprefix}_hwsex
    xcount=$( awk '$1 == 23' ${prunehqprefix}_hwsex.bim | wc -l )
    if (( xcount > minvarcount )) ; then
      plink --bfile ${prunehqprefix}_isex --extract <( cut -f2 ${prunehqprefix}_hwsex.bim ) \
        --impute-sex --make-bed --out ${prunehqprefix}_isex
    fi
    usex="--update-sex ${prunehqprefix}_isex.fam 3"
    awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' \
      ${prunehqprefix}_isex.fam > ${prunehqprefix}_isex.nosex
    if [[ -s "${prunehqprefix}_isex.nosex" ]] ; then
      nosex="--remove ${prunehqprefix}_isex.nosex"
    fi
  fi
  tmpbiofile=$( mktemp ${mybiofile}.tmpXXXX )
  sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;' ${prunehqprefix}_isex.sexcheck \
  | awk -F $'\t' -v uid=${uid} '{
    OFS="\t"
    if ( NR>1 ) uid=$1"_"$2
    printf( "%s", uid )
    for ( k=3; k<=NF; k++ )
      printf( "\t%s", $k )
    printf( "\n" )
  }' | sort -t $'\t' -u -k 1,1 | join -t $'\t' -a1 -e '-' ${mybiofile} - > ${tmpbiofile}
  mv ${tmpbiofile} ${mybiofile}
fi
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

