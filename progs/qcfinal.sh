#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r  tmpprefix=${opt_outprefix}_tmp
declare -r  debuglogfn=${tmpprefix}_debug.log
declare -ra batchfamfiles=( $( ls ${opt_batchoutprefix}*.fam ) )

declare -r cfg_hweneglogp_ctrl=4
declare -r cfg_minindcount=100
declare -r cfg_phenotypes='/cluster/projects/p33/users/franbe/norment_2018/test/phenotypes.txt'
declare -r cfg_samplemiss=0.05

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping final QC step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
  exit 1
fi

# input: individual qc'd plink set
# output: plink genotype set for variant passing final QC criteria:
#         - more stringent control Hardy-Weinberg equilibrium
#         - no control batch effects


# if there are any annotated controls in the original file use those
cut -f 1,2,6 ${opt_inprefix}.fam | sed -r 's/[ \t]+/\t/g' | sort -u > ${tmpprefix}_ctrl.txt

declare Nctrl=0
if [ "${cfg_phenotypes}" != "" -a -s "${cfg_phenotypes}" ] ; then
  Nctrl=$( awk -F $'\t' 'tolower($3) == "control"' ${cfg_phenotypes} | wc -l )
fi
declare plinkflag=''
if [ $Nctrl -ge $cfg_minindcount ] ; then
  # overwrite control set with the provided custom list
  awk -F $'\t' 'tolower($3) == "control"' ${cfg_phenotypes} | sort -u > ${tmpprefix}_ctrl.txt
  # redefine phenotypes in the input plink set
  plink --bfile ${opt_inprefix} \
        --make-pheno ${cfg_phenotypes} affected \
        --make-bed \
        --out ${tmpprefix}_pheno \
        >> ${debuglogfn}
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_pheno.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_pheno.fam
  # enforce stricter non-sex chromosomes Hardy-Weinberg equilibrium on controls
  plink --bfile ${tmpprefix}_pheno ${nosex_flag} \
        --not-chr 23,24 \
        --hwe 1.E-${cfg_hweneglogp_ctrl} \
        --make-just-bim \
        --out ${tmpprefix}_ctrlhwe_nonsex \
        >> ${debuglogfn}
  nosex_flag=''
  # get set of individuals missing sex information for exclusion from later check
  awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' ${opt_hqprefix}.fam \
    > ${opt_hqprefix}.nosex
  if [ -s "${opt_hqprefix}.nosex" ] ; then
    nosex_flag="--remove ${opt_hqprefix}.nosex"
  fi
  sex_hweneglogp_ctrl=$(( cfg_hweneglogp_ctrl/2 ))
  if [ "${sex_hweneglogp_ctrl}" -lt 12 ] ; then
    sex_hweneglogp_ctrl=12
  fi
  # enforce stricter sex chromosomes Hardy-Weinberg equilibrium on controls
  plink --bfile ${tmpprefix}_pheno ${nosex_flag} \
        --chr 23,24 \
        --hwe 1.E-${sex_hweneglogp_ctrl} \
        --make-just-bim \
        --out ${tmpprefix}_ctrlhwe_sex \
        >> ${debuglogfn}
  [ -s "${tmpprefix}_ctrlhwe_nonsex.bim" -o -s "${tmpprefix}_ctrlhwe_sex.bim" ] || {
    printf "error: no variants left after HWE.\n" >&2;
    exit 1;
  }
  # list the variants passing stricter Hardy-Weiberg equilibrium tests on controls
  cut -f 2 ${tmpprefix}_ctrlhwe_*sex.bim | sort -u > ${tmpprefix}_ctrlhwe.mrk
  plinkflag="--extract ${tmpprefix}_ctrlhwe.mrk"
fi

# make a new plink set with eventual filters
plink --bfile ${opt_inprefix} ${plinkflag} \
  --make-bed \
  --out ${tmpprefix}_ctrlhwe \
  >> ${debuglogfn}
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_ctrlhwe.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_ctrlhwe.fam

# make a copy of the files with output suffix
cp ${tmpprefix}_ctrlhwe.bed ${tmpprefix}_out.bed
cp ${tmpprefix}_ctrlhwe.bim ${tmpprefix}_out.bim
cp ${tmpprefix}_ctrlhwe.fam ${tmpprefix}_out.fam

if [ ${#batchfamfiles[*]} -gt 1 ] ; then
  > ${tmpprefix}.exclude
  for i in ${!batchfamfiles[@]} ; do
    echo "assessing batch effects for '${batchfamfiles[$i]}'.."
    plink --bfile ${opt_inprefix} \
          --keep ${tmpprefix}_ctrl.txt \
          --make-pheno ${batchfamfiles[$i]} '*' \
          --model \
          --out ${tmpprefix}_plink \
          >> ${debuglogfn}
    if [ -s "${tmpprefix}_plink.model" ] ; then
      sed -i -r 's/^[ \t]+//g; s/[ \t]+/\t/g;' ${tmpprefix}_plink.model
      for atest in $( cut -f 5 ${tmpprefix}_plink.model | tail -n +2 | sort -u ) ; do
        echo "summarizing ${atest} tests.."
        awk -F $'\t' -v atest=${atest} 'NR > 1 && $5 == atest' ${tmpprefix}_plink.model \
          | cut -f 2,10 \
          | ${BASEDIR}/progs/fdr.Rscript \
          | awk -F $'\t' '$2 < 0.9' \
          | cut -f 1 \
          >> ${tmpprefix}.exclude
      done
    fi
  done
  sort -u ${tmpprefix}.exclude > ${tmpprefix}.exclude.sort
  mv ${tmpprefix}.exclude.sort ${tmpprefix}.exclude
  plink --bfile ${tmpprefix}_ctrlhwe \
        --exclude ${tmpprefix}.exclude \
        --make-bed \
        --out ${tmpprefix}_out \
        >> ${debuglogfn}
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
fi

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
# purge fam file ids of any unwanted characters
cp ${tmpprefix}_out.fam ${opt_outprefix}.fam.org
awk '{
  OFS="\t"
  for ( k = 1; k < 5; k++ )
    gsub( "[/+-]", "_", $k )
  print
}' ${tmpprefix}_out.fam > ${opt_outprefix}.fam

