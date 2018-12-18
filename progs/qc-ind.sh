#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_samplemiss=$( cfgvar_get samplemiss )

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping individual QC step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

# input: variant qc'd plink set
# output: plink genotype set for individuals passing QC criteria:
#         - good coverage


tmp_samplemiss=${cfg_samplemiss}
N=$( wc -l ${opt_inprefix}.bim | cut -d ' ' -f 1 )
if [ $N -lt ${cfg_minvarcount} ] ; then tmp_samplemiss=0.1 ; fi
${plinkexec} --bfile ${opt_inprefix} \
             --set-hh-missing \
             --set-me-missing \
             --mind ${tmp_samplemiss} \
             --make-bed \
             --out ${tmpprefix}_out \
             2>&1 >> ${debuglogfn} \
             | tee -a ${debuglogfn}
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
{
  awk -F $'\t' '{ print( $1"_"$2 ); }' ${tmpprefix}_out.fam \
    | sort -u \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' '{
      OFS="\t"
      if ( NR == 1 ) print( $0, "COV" )
      else {
        if ( $(NF-2) == "PROBLEM" ) print( $0, 0 )
        else print( $0, "NO" )
      }
    }'
  awk -F $'\t' '{ print( $1"_"$2 ); }' ${tmpprefix}_out.fam \
    | sort -u \
    | join -t $'\t' ${opt_biofile} - \
    | awk -F $'\t' '{
      OFS="\t"
      print( $0, "YES" )
    }'
} | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.bio
cp ${tmpprefix}.bio ${opt_biofile}

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
mv ${tmpprefix}_out.log ${opt_outprefix}.log

rm ${tmpprefix}*

