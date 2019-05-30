#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare -r  tmpprefix="${opt_outprefix}_tmp"
declare -ra batchfamfiles=( $( ls "${opt_batchoutprefix}"*.fam ) )

declare -r cfg_bfdr=$( cfgvar_get bfdr )
declare -r cfg_hweneglogp_ctrl=$( cfgvar_get hweneglogp_ctrl )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_phenotypes=$( cfgvar_get phenotypes )

#-------------------------------------------------------------------------------

# input: variant qc'd plink set
# output: plink genotype set for variant passing final QC criteria:

printf "\
  * Extract variants with:
    * more stringent control Hardy-Weinberg equilibrium
    * eventual batch effects (FDR=${cfg_bfdr}) (if more than a single batch present)
\n" | printlog 1

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "> '%s' found. skipping final QC..\n" "${opt_outprefix}.bed"
  exit 0
fi

if ls "${tmpprefix}"* > /dev/null 2>&1; then
  printf "> temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

declare keepfile=''
declare keepflag=''
# set keep flag if a list of unrelated individuals exists
if [ -f "${opt_outprefixbase}/.i/qc/e_indqc.ids" ] ; then
  keepfile="${opt_outprefixbase}/.i/qc/e_indqc.ids"
  keepflag="--keep ${keepfile}"
fi

cp "${opt_inprefix}.bed" "${tmpprefix}_out.bed"
cp "${opt_inprefix}.bim" "${tmpprefix}_out.bim"
cp "${opt_inprefix}.fam" "${tmpprefix}_out.fam"

if [ "${cfg_phenotypes}" != "" -a -s "${cfg_phenotypes}" ] ; then
  # if a phenotype file was specified write a control list
  printf "> extracting control list from '${cfg_phenotypes}'..\n"
  awk -F $'\t' '$3 == 1' ${cfg_phenotypes} \
    | sort -u \
    | extract_sample_ids "${keepfile}" \
    > "${tmpprefix}_ctrl.txt"
  declare -r Nctrl=$( cat "${tmpprefix}_ctrl.txt" | wc -l )
  # and redefine phenotypes in the input plink set
  ${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" \
               --make-pheno "${cfg_phenotypes}" 2 \
               --make-bed \
               --out "${tmpprefix}_pheno" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
elif [ $( grep -c '1$' "${opt_inprefix}.fam" ) -ge $cfg_minindcount ] ; then
  # else, if there are enough annotated controls in the original file use those
  printf "> extracting control list from '${opt_inprefix}.fam'..\n"
  awk -F $'\t' '$6 == 1' "${opt_inprefix}.fam" \
    | cut -f 1,2,6 \
    | sort -u \
    | extract_sample_ids "${keepfile}" \
    > "${tmpprefix}_ctrl.txt"
  declare -r Nctrl=$( cat "${tmpprefix}_ctrl.txt" | wc -l )
  # rename temporary plink set
  rename _out _pheno "${tmpprefix}_out".*
else
  # else, use the whole list but leave Nctrl=0 to suppress control-HWE tests
  printf "> no controls available. using everyone..\n"
  cut -f 1,2,6 "${opt_inprefix}.fam" \
    | extract_sample_ids "${keepfile}" \
    | sort -u \
    > "${tmpprefix}_ctrl.txt" 
  declare -r Nctrl=0
fi
printf "> $Nctrl actual controls found.\n"
# if Nctrl is not zero a plink *_pheno set should exist
if [ $Nctrl -ge $cfg_minindcount ] ; then
  declare plinkflag=''
  # enforce stricter non-sex chromosomes Hardy-Weinberg equilibrium on controls
  printf "> testing non-sex chromosomes control Hardy-Weinberg equilibrium..\n"
  echo "  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_pheno ${keepflag}
               --not-chr 23,24
               --hwe 1.E-${cfg_hweneglogp_ctrl} midp
               --make-just-bim
               --out ${tmpprefix}_ctrlhwe_nonsex" | printlog 2
  ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_pheno" ${keepflag} \
               --not-chr 23,24 \
               --hwe 1.E-${cfg_hweneglogp_ctrl} midp \
               --make-just-bim \
               --out "${tmpprefix}_ctrlhwe_nonsex" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  declare nosexflag=''
  # get set of individuals missing sex information for exclusion from later check
  awk '{ OFS="\t"; if ( NR > 1 && $5 == 0 ) print( $1, $2 ); }' "${opt_hqprefix}.fam" \
    > "${opt_hqprefix}.nosex"
  if [ -s "${opt_hqprefix}.nosex" ] ; then
    nosexflag="--remove ${opt_hqprefix}.nosex"
  fi
  sex_hweneglogp_ctrl=$(( cfg_hweneglogp_ctrl*2 ))
  if [ "${sex_hweneglogp_ctrl}" -gt 12 ] ; then
    sex_hweneglogp_ctrl=12
  fi
  if [ $( get_xvar_count "${tmpprefix}_pheno.bim" ) -ge ${cfg_minvarcount} ] ; then
    # enforce stricter sex chromosomes Hardy-Weinberg equilibrium on controls
    printf "> testing sex chromosomes control Hardy-Weinberg equilibrium..\n"
    echo "  ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_pheno" ${keepflag} ${nosexflag}
                 --chr 23,24
                 --hwe 1.E-${sex_hweneglogp_ctrl} midp
                 --make-just-bim
                 --out ${tmpprefix}_ctrlhwe_sex" | printlog 2
    ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_pheno" ${keepflag} ${nosexflag} \
                 --chr 23,24 \
                 --hwe 1.E-${sex_hweneglogp_ctrl} midp \
                 --make-just-bim \
                 --out "${tmpprefix}_ctrlhwe_sex" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
  fi
  unset nosexflag
  [ -s "${tmpprefix}_ctrlhwe_nonsex.bim" -o -s "${tmpprefix}_ctrlhwe_sex.bim" ] || {
    printf "error: no variants left after HWE.\n" >&2;
    exit 1;
  }
  # list the variants passing stricter Hardy-Weiberg equilibrium tests on controls
  cut -f 2 "${tmpprefix}_ctrlhwe_"*sex.bim | sort -u > "${tmpprefix}_ctrlhwe.mrk"
  plinkflag="--extract ${tmpprefix}_ctrlhwe.mrk"
  # make a new plink set with eventual filter
  ${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" ${plinkflag} \
               --make-bed \
               --out "${tmpprefix}_ctrlhwe" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # make a copy of the files with output suffix
  cp "${tmpprefix}_ctrlhwe.bed" "${tmpprefix}_out.bed"
  cp "${tmpprefix}_ctrlhwe.bim" "${tmpprefix}_out.bim"
  cp "${tmpprefix}_ctrlhwe.fam" "${tmpprefix}_out.fam"
fi
printf "> ${#batchfamfiles[*]} batches found.\n"
if [ ${#batchfamfiles[*]} -gt 1 ] ; then
  > "${tmpprefix}.exclude"
  for i in ${!batchfamfiles[@]} ; do
    printf "> assessing batch effects for '${batchfamfiles[$i]}'..\n"
    ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_out" \
                 --allow-no-sex \
                 --keep "${tmpprefix}_ctrl.txt" \
                 --make-pheno ${batchfamfiles[$i]} '*' \
                 --model \
                 --out "${tmpprefix}_plink" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    if [ -s "${tmpprefix}_plink.model" ] ; then
      sed -i -r 's/^[ \t]+//g; s/[ \t]+/\t/g;' "${tmpprefix}_plink.model"
#       awk '$10 == "NA"' "${tmpprefix}_plink.model" >> "${tmpprefix}.exclude"
      for atest in $( cut -f 5 "${tmpprefix}_plink.model" | tail -n +2 | sort -u ) ; do
        printf "> summarizing ${atest} tests..\n"
        n=$( awk -F $'\t' -v a=${atest} 'NR > 1 && $5 == a && $10 != "NA"' \
          "${tmpprefix}_plink.model" | wc -l )
        tail -n +2 "${tmpprefix}_plink.model" | sort -t $'\t' -k 10,10gr \
          | awk -F $'\t' -f "${BASEDIR}/lib/awk/stats.awk" \
            -v a=${atest} -v bfdr=${cfg_bfdr} -v n=${n} \
            --source 'BEGIN{
                cnt = 0
                BHfdr = 1
              } {
                if ( $5 == a && $10 != "NA" ) {
                  BHfdr = pmin(BHfdr, $10*(n/(n - cnt + 1)));
                  if ( BHfdr < bfdr ) print $2
                  cnt++
                }
              }' >> "${tmpprefix}.exclude"
      done
    fi
  done
  sort -u "${tmpprefix}.exclude" > "${tmpprefix}.exclude.sort"
  mv "${tmpprefix}.exclude.sort" "${tmpprefix}.exclude"
  ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_out" \
               --exclude "${tmpprefix}.exclude" \
               --make-bed \
               --out "${tmpprefix}_nbe" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # make a copy of the files with output suffix
  cp "${tmpprefix}_nbe.bed" "${tmpprefix}_out.bed"
  cp "${tmpprefix}_nbe.bim" "${tmpprefix}_out.bim"
  cp "${tmpprefix}_nbe.fam" "${tmpprefix}_out.fam"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"

rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out".*

rm "${tmpprefix}"*

