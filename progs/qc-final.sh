#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r  tmpprefix=${opt_outprefix}_tmp
declare -r  debuglogfn=${tmpprefix}_debug.log
declare -ra batchfamfiles=( $( ls ${opt_batchoutprefix}*.fam ) )

declare -r cfg_bfdr=$( cfgvar_get bfdr )
declare -r cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_phenotypes=$( cfgvar_get phenotypes )

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping final QC step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

# input: individual qc'd plink set
# output: plink genotype set for variant passing final QC criteria:
#         - more stringent control Hardy-Weinberg equilibrium
#         - no control batch effects

echo 'performing final qc..'

declare keepfile=''
declare keepflag=''
# set keep flag if a list of unrelated individuals exists
if [ -f "${opt_outprefixbase}.ids" ] ; then
  keepfile=${opt_outprefixbase}.ids
  keepflag="--keep ${keepfile}"
fi

cp ${opt_inprefix}.bed ${tmpprefix}_out.bed
cp ${opt_inprefix}.bim ${tmpprefix}_out.bim
cp ${opt_inprefix}.fam ${tmpprefix}_out.fam

if [ "${cfg_phenotypes}" != "" -a -s "${cfg_phenotypes}" ] ; then
  # if a phenotype file was specified write a control list
  echo "extracting control list from '${cfg_phenotypes}'.."
  awk -F $'\t' '$3 == 1' ${cfg_phenotypes} \
    | sort -u \
    | extract_sample_ids ${keepfile} \
    > ${tmpprefix}_ctrl.txt
  declare -r Nctrl=$( cat ${tmpprefix}_ctrl.txt | wc -l )
  # and redefine phenotypes in the input plink set
  ${plinkexec} --bfile ${opt_inprefix} \
               --make-pheno ${cfg_phenotypes} 2 \
               --make-bed \
               --out ${tmpprefix}_pheno \
               2>&1 >> ${debuglogfn} \
               | tee -a ${debuglogfn}
elif [ $( grep -c '1$' ${opt_inprefix}.fam ) -ge $cfg_minindcount ] ; then
  # else, if there are enough annotated controls in the original file use those
  echo "extracting control list from '${opt_inprefix}.fam'.."
  awk -F $'\t' '$6 == 1' ${opt_inprefix}.fam \
    | cut -f 1,2,6 \
    | sort -u \
    | extract_sample_ids ${keepfile} \
    > ${tmpprefix}_ctrl.txt
  declare -r Nctrl=$( cat ${tmpprefix}_ctrl.txt | wc -l )
  # rename temporary plink set
  rename _out _pheno ${tmpprefix}_out.*
else
  # else, use the whole list but leave Nctrl=0 to suppress control-HWE tests
  echo "no controls available. using everyone.."
  cut -f 1,2,6 ${opt_inprefix}.fam \
    | extract_sample_ids ${keepfile} \
    | sort -u \
    > ${tmpprefix}_ctrl.txt 
  declare -r Nctrl=0
fi
echo "$Nctrl actual controls found."
# if Nctrl is not zero a plink *_pheno set should exist
if [ $Nctrl -ge $cfg_minindcount ] ; then
  declare plinkflag=''
  # enforce stricter non-sex chromosomes Hardy-Weinberg equilibrium on controls
  echo "testing non-sex chromosomes control Hardy-Weinberg equilibrium.."
  ${plinkexec} --bfile ${tmpprefix}_pheno ${keepflag} \
               --not-chr 23,24 \
               --hwe 1.E-${cfg_hweneglogp_ctrl} midp \
               --make-just-bim \
               --out ${tmpprefix}_ctrlhwe_nonsex \
               2>&1 >> ${debuglogfn} \
               | tee -a ${debuglogfn}
  declare nosexflag=''
  # get set of individuals missing sex information for exclusion from later check
  awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' ${opt_hqprefix}.fam \
    > ${opt_hqprefix}.nosex
  if [ -s "${opt_hqprefix}.nosex" ] ; then
    nosexflag="--remove ${opt_hqprefix}.nosex"
  fi
  sex_hweneglogp_ctrl=$(( cfg_hweneglogp_ctrl*2 ))
  if [ "${sex_hweneglogp_ctrl}" -gt 12 ] ; then
    sex_hweneglogp_ctrl=12
  fi
  if [ $( get_xvar_count ${tmpprefix}_pheno.bim ) -ge ${cfg_minvarcount} ] ; then
    # enforce stricter sex chromosomes Hardy-Weinberg equilibrium on controls
    echo "testing sex chromosomes control Hardy-Weinberg equilibrium.."
    ${plinkexec} --bfile ${tmpprefix}_pheno ${keepflag} ${nosexflag} \
                 --chr 23,24 \
                 --hwe 1.E-${sex_hweneglogp_ctrl} midp \
                 --make-just-bim \
                 --out ${tmpprefix}_ctrlhwe_sex \
                 2>&1 >> ${debuglogfn} \
                 | tee -a ${debuglogfn}
  fi
  unset nosexflag
  [ -s "${tmpprefix}_ctrlhwe_nonsex.bim" -o -s "${tmpprefix}_ctrlhwe_sex.bim" ] || {
    printf "error: no variants left after HWE.\n" >&2;
    exit 1;
  }
  # list the variants passing stricter Hardy-Weiberg equilibrium tests on controls
  cut -f 2 ${tmpprefix}_ctrlhwe_*sex.bim | sort -u > ${tmpprefix}_ctrlhwe.mrk
  plinkflag="--extract ${tmpprefix}_ctrlhwe.mrk"
  # make a new plink set with eventual filter
  ${plinkexec} --bfile ${opt_inprefix} ${plinkflag} \
               --make-bed \
               --out ${tmpprefix}_ctrlhwe \
               2>&1 >> ${debuglogfn} \
               | tee -a ${debuglogfn}
  # make a copy of the files with output suffix
  cp ${tmpprefix}_ctrlhwe.bed ${tmpprefix}_out.bed
  cp ${tmpprefix}_ctrlhwe.bim ${tmpprefix}_out.bim
  cp ${tmpprefix}_ctrlhwe.fam ${tmpprefix}_out.fam
fi

echo "${#batchfamfiles[*]} batches found."
if [ ${#batchfamfiles[*]} -gt 1 ] ; then
  > ${tmpprefix}.exclude
  for i in ${!batchfamfiles[@]} ; do
    echo "assessing batch effects for '${batchfamfiles[$i]}'.."
    ${plinkexec} --bfile ${tmpprefix}_out \
                 --allow-no-sex \
                 --keep ${tmpprefix}_ctrl.txt \
                 --make-pheno ${batchfamfiles[$i]} '*' \
                 --model \
                 --out ${tmpprefix}_plink \
                 2>&1 >> ${debuglogfn} \
                 | tee -a ${debuglogfn}
    if [ -s "${tmpprefix}_plink.model" ] ; then
      sed -i -r 's/^[ \t]+//g; s/[ \t]+/\t/g;' ${tmpprefix}_plink.model
      for atest in $( cut -f 5 ${tmpprefix}_plink.model | tail -n +2 | sort -u ) ; do
        echo "summarizing ${atest} tests.."
        awk -F $'\t' -v atest=${atest} 'NR > 1 && $5 == atest' ${tmpprefix}_plink.model \
          | cut -f 2,10 \
          | ${BASEDIR}/progs/fdr.Rscript \
          | awk -v bfdr=${cfg_bfdr} -F $'\t' '$2 < bfdr' \
          | cut -f 1 \
          >> ${tmpprefix}.exclude
      done
    fi
  done
  sort -u ${tmpprefix}.exclude > ${tmpprefix}.exclude.sort
  mv ${tmpprefix}.exclude.sort ${tmpprefix}.exclude
  ${plinkexec} --bfile ${tmpprefix}_out \
               --exclude ${tmpprefix}.exclude \
               --make-bed \
               --out ${tmpprefix}_nbe \
               2>&1 >> ${debuglogfn} \
               | tee -a ${debuglogfn}
  # make a copy of the files with output suffix
  cp ${tmpprefix}_nbe.bed ${tmpprefix}_out.bed
  cp ${tmpprefix}_nbe.bim ${tmpprefix}_out.bim
  cp ${tmpprefix}_nbe.fam ${tmpprefix}_out.fam
fi

sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
# purge fam file ids of any unwanted characters
cp ${tmpprefix}_out.fam ${opt_outprefix}.fam.org
awk '{
  OFS="\t"
  for ( k = 1; k < 5; k++ )
    gsub( "[][)(}{/\\,.;:|!?@#$%^&*~=_><+-]+", "_", $k )
  print
}' ${opt_outprefix}.fam.org > ${tmpprefix}_out.fam
Nold=$( sort -u -k 1,2 ${opt_outprefix}.fam.org | wc -l )
Nnew=$( sort -u -k 1,2 ${tmpprefix}_out.fam | wc -l )
if [ $Nold -eq $Nnew ] ; then
  mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
else
  echo '====> warning: could not purge IDs due to conflicts.'
  mv ${opt_outprefix}.fam.org ${opt_outprefix}.fam
fi

rm ${tmpprefix}*

