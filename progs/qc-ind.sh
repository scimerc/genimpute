#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare -r tmpprefix="${opt_outprefix}_tmp"

declare -r cfg_hvm=$( cfgvar_get hvm )
declare -r cfg_hvmmax=$( cfgvar_get hvmmax )
declare -r cfg_hvmmin=$( cfgvar_get hvmmin )
declare -r cfg_hvmsdn=$( cfgvar_get hvmsdn )
declare -r cfg_pihatrel=$( cfgvar_get pihatrel )
declare -r cfg_minindcount=$( cfgvar_get minindcount )
declare -r cfg_minvarcount=$( cfgvar_get minvarcount )
declare -r cfg_samplemiss=$( cfgvar_get samplemiss )
declare -r cfg_varmiss=$( cfgvar_get varmiss )
declare -r cfg_uid=$( cfgvar_get uid )

declare -r kingexec="${BASEDIR}/lib/3rd/king"

#-------------------------------------------------------------------------------

break_dysfunc_fam() {
  # in: plink fam file
  # out: plink update-parents file
  cut -f 1 "$@" | sort | uniq -c | tabulate \
    | awk -F $'\t' '{ OFS="\t"; if ( $1 > 1 ) print( $2 ); }' \
    | join -t $'\t' - <( sort -k 1,1 "$@" ) \
    | awk '{
        OFS="\t"
        dad[$1,$2] = $1 SUBSEP $3
        mom[$1,$2] = $1 SUBSEP $4
        sex[$1,$2] = $5
      } END{
        for ( sid in sex ) {
          parsex = sex[dad[sid]] == 0 || sex[mom[sid]] == 0
          dysfam = parsex || sex[dad[sid]] == sex[mom[sid]]
          split( sid, ovid, SUBSEP )
          split( mom[sid], mvid, SUBSEP )
          split( dad[sid], dvid, SUBSEP )
          if ( dysfam ) print( ovid[1], ovid[2], 0, 0 )
          else print( ovid[1], ovid[2], dvid[2], mvid[2] )
        }
      }'
}

#-------------------------------------------------------------------------------

# input: merged plink set and hq plink set
# output: clean plink set (with imputed sex from hq set and no mixups)

echo -e "==== Individual QC ====\n" | printlog 1

printf "\
  * Exclude possibly contaminated (mean+%dsd heterozygosity) and low-coverage (%.0f%%) individuals
  * Erase parent information for any dysfunctional families
  * Annotate relatedness information
\n" ${cfg_hvmsdn} "$( echo "(1 - ${cfg_hvmmax})*100" | bc )"  | printlog 1

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "> '%s' found. skipping individual QC..\n" "${opt_outprefix}.bed"
  exit 0
fi

if ls "${tmpprefix}"* > /dev/null 2>&1; then
  printf "> temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

#-------------------------------------------------------------------------------

declare keepflag=''
# find potential mixups?
if [ ${cfg_hvm} -eq 1 ] ; then
  printf "> computing individual heterozygosity and missing rates..\n"
  ${plinkexec} --allow-extra-chr \
               --bfile "${opt_hqprefix}i" \
               --het --missing \
               --out "${tmpprefix}_sq" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  printf "> removing mixups and low-coverage individuals..\n"
  join --header -t $'\t' \
    <( attach_uids -h "${tmpprefix}_sq.imiss" | sort -t $'\t' -u -k 1,1 ) \
    <( attach_uids -h "${tmpprefix}_sq.het" | cut -f 1,4- | sort -t $'\t' -u -k 1,1 ) \
    | awk -F $'\t' \
      -f "${BASEDIR}/lib/awk/stats.awk" \
      -v maxmis=${cfg_hvmmax} \
      -v minmis=${cfg_hvmmin} \
      -v sdn=${cfg_hvmsdn} \
      --source '{
        if ( NR > 1 ) {
          fid[$1] = $2
          iid[$1] = $3
          het[$1] = ( $10 - $8 ) / $10
          mis[$1] = $7
        }
      } END{
        OFS = "\t"
        meanhet = smean( het )
        sdhet = ssd( het )
        meanmis = smean( mis )
        sdmis = ssd( mis )
        misthreshmax = maxmis
        misthresh = pmin( meanmis + sdn*sdmis, minmis )
        for ( uid in iid ) {
          # keep samples with low heterozygosity or too low coverage (within limits) to tell
          goodhet = mis[uid] <= misthresh && het[uid] <= meanhet + sdn*sdhet
          highmis = mis[uid] >= misthresh && mis[uid] <= misthreshmax
          if ( goodhet || highmis ) print( fid[uid], iid[uid] )
        }
      }' | sort -u \
      > "${tmpprefix}_out.clean.id"
  [ -s "${tmpprefix}_out.clean.id" ] || {
    printf "> error: file '%s' empty or not found.\n" "${tmpprefix}_out.clean.id" >&2;
    exit 1;
  }
  # write plink flag for non-mixup info later
  keepflag="--keep ${tmpprefix}_out.clean.id"
else
  # write a dummy clean.id file including everyone from the original file
  cut -f 1,2 "${opt_inprefix}.fam" | sort -u > "${tmpprefix}_out.clean.id"
fi
# extract clean, high coverage individuals
printf "> extracting high coverage individuals..\n"
break_dysfunc_fam "${opt_hqprefix}.fam" > "${tmpprefix}_out_updateparents.txt"
tmp_samplemiss=${cfg_samplemiss}
N=$( wc -l "${opt_hqprefix}.bim" | cut -d ' ' -f 1 )
if [ $N -lt ${cfg_minvarcount} ] ; then tmp_samplemiss=0.1 ; fi
echo "  ${plinkexec} --allow-extra-chr \
             --bfile ${opt_hqprefix} ${keepflag}
             --mind ${tmp_samplemiss}
             --update-parents ${tmpprefix}_out_updateparents.txt
             --make-bed
             --out ${tmpprefix}_hc" | printlog 2
${plinkexec} --allow-extra-chr \
             --bfile "${opt_hqprefix}" ${keepflag} \
             --mind ${tmp_samplemiss} \
             --update-parents "${tmpprefix}_out_updateparents.txt" \
             --make-bed \
             --out "${tmpprefix}_hc" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
tmp_samplemiss=${cfg_samplemiss}
N=$( wc -l "${opt_hqprefix}i.bim" | cut -d ' ' -f 1 )
if [ $N -lt ${cfg_minvarcount} ] ; then tmp_samplemiss=0.1 ; fi
echo "  ${plinkexec} --allow-extra-chr \
             --bfile ${opt_hqprefix}i ${keepflag}
             --mind ${tmp_samplemiss}
             --update-parents ${tmpprefix}_out_updateparents.txt
             --make-just-fam
             --out ${tmpprefix}_hci" | printlog 2
${plinkexec} --allow-extra-chr \
             --bfile "${opt_hqprefix}i" ${keepflag} \
             --mind ${tmp_samplemiss} \
             --update-parents "${tmpprefix}_out_updateparents.txt" \
             --make-just-fam \
             --out "${tmpprefix}_hci" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# reconstruct pedigrees using king
printf "> reconstructing eventual families (second degree)..\n"
${kingexec}  -b "${tmpprefix}_hc.bed" \
             --build --degree 2 --prefix "${tmpprefix}_out_" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# write expected king output files if these are missing or empty
[ -s "${tmpprefix}_out_updateids.txt" ] || awk '{ print( $1, $2, $1, $2 ); }' \
  "${tmpprefix}_hc.fam" > "${tmpprefix}_out_updateids.txt"
[ -s "${tmpprefix}_out_updateparents.txt" ] || awk '{ print( $1, $2, $3, $4 ); }' \
  "${tmpprefix}_hc.fam" > "${tmpprefix}_out_updateparents.txt"
${plinkexec} --allow-extra-chr \
             --fam "${tmpprefix}_hc.fam" \
             --update-ids "${tmpprefix}_out_updateids.txt" \
             --make-just-fam \
             --out "${tmpprefix}_kingids" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
${plinkexec} --allow-extra-chr \
             --fam "${tmpprefix}_kingids.fam" \
             --update-parents "${tmpprefix}_out_updateparents.txt" \
             --make-just-fam \
             --out "${tmpprefix}_kingpeds" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_king"*.fam
# reannotate eventual dysfunctional family information
printf "> reannotating eventual dysfunctional families..\n"
break_dysfunc_fam "${tmpprefix}_kingpeds.fam" > "${tmpprefix}_out_updateparents.txt"
# identify related individuals
printf "> identifying related individuals..\n"
${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}i" \
             --keep "${tmpprefix}_hci.fam" \
             --genome gz \
             --out "${tmpprefix}_sq" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}i" \
             --keep "${tmpprefix}_hci.fam" \
             --cluster \
             --read-genome "${tmpprefix}_sq.genome.gz" \
             --rel-cutoff ${cfg_pihatrel} \
             --out "${tmpprefix}_sq" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# rename list of clean, unrelated individuals for later use
mv "${tmpprefix}_sq.rel.id" "${opt_outprefixbase}/.i/qc/e_indqc.ids"
# remove mixups, update sex and set heterozygous haploid genotypes to missing in input set
${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" ${keepflag} \
             --update-sex "${opt_hqprefix}.fam" 3 \
             --set-hh-missing \
             --make-bed \
             --out "${tmpprefix}_out" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"
unset keepflag

extract_related_lists_from_grm_file() {
  local -r infile="$1"
  zcat -f "${infile}" | tabulate \
    | awk -F $'\t' -v uid=${cfg_uid} -v pihat=${cfg_pihatrel} \
      -f "${BASEDIR}"/lib/awk/idclean.awk --source '{
        maxcnt=11111
        uid0=idclean( $1"_"$2 )
        uid1=idclean( $3"_"$4 )
        if ( NR>1 && $10>=pihat ) {
          if ( uid0 in relarr && cntarr[uid0] < maxcnt ) {
            relarr[uid0] = relarr[uid0]","uid1"("$10")"
            cntarr[uid0]++
          }
          else {
            relarr[uid0] = uid1"("$10")"
            cntarr[uid0] = 1
          }
          if ( uid1 in relarr && cntarr[uid1] < maxcnt ) {
            relarr[uid1] = relarr[uid1]","uid0"("$10")"
            cntarr[uid1]++
          }
          else {
            relarr[uid1] = uid0"("$10")"
            cntarr[uid1] = 1
          }
        }
      } END{
        OFS="\t"
        printf( "%s\tRELSHIP\n", uid )
        for ( uid in relarr ) print( uid, relarr[uid] )
      }' \
    | sort -t $'\t' -u -k 1,1
}

# update biography file with sex information
printf "> updating biography file with sex information..\n"
{
  {
    printf "%s\t" ${cfg_uid}
    cat "${opt_hqprefix}"*.sexcheck > "${tmpprefix}.sexcheck"
    head -n 1 "${tmpprefix}.sexcheck" | tabulate | cut -f 3-
  } | join -t $'\t'     "${opt_biofile}" - \
    | tee "${tmpprefix}.bio.head"
  # count number of fields in the merged file
  TNF=$( cat "${tmpprefix}.bio.head" | wc -w )
  # merge information for existing individuals
  attach_uids "${opt_hqprefix}"*.sexcheck \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t'     "${opt_biofile}" - \
  # add non-existing individuals and pad the extra fields with NAs
  attach_uids "${opt_hqprefix}"*.sexcheck \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | awk -F $'\t' '{
    OFS="\t"
    if ( NR==1 ) {
      for ( k=1; k<=NF; k++ )
        if ( $k ~ "^F$" ) $k = "F_SEX"
    }
    print
  }' | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.0.bio"
mv "${tmpprefix}.0.bio" "${opt_biofile}"

# update biography file with heterozygosity information
printf "> updating biography file with heterozygosity information..\n"
{
  {
    printf "%s\t" ${cfg_uid}
    head -n 1 "${tmpprefix}_sq.het" | tabulate | cut -f 3-
  } | join -t $'\t'     "${opt_biofile}" - \
    | tee "${tmpprefix}.bio.head"
  # count number of fields in the merged file
  TNF=$( cat "${tmpprefix}.bio.head" | wc -w )
  # merge information for existing individuals
  attach_uids "${tmpprefix}_sq.het" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t'     "${opt_biofile}" - \
  # add non-existing individuals and pad the extra fields with NAs
  attach_uids "${tmpprefix}_sq.het" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | awk -F $'\t' '{
    OFS="\t"
    if ( NR==1 ) {
      for ( k=1; k<=NF; k++ )
        if ( $k ~ "^F$" ) $k = "F_HET"
    }
    print
  }' | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.1.bio"
mv "${tmpprefix}.1.bio" "${opt_biofile}"

# update biography file with potential mixup information
printf "> updating biography file with potential contaminations..\n"
{
  attach_uids "${tmpprefix}_out.clean.id" \
    | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' -v hvm=${cfg_hvm} '{
      OFS="\t"
      hvmtag="PROBLEM"
      if ( hvm != 1 ) hvmtag="__NA__"
      if ( NR == 1 ) print( $0, "MISMIX" )
      else print( $0, hvmtag )
    }'
  attach_uids "${tmpprefix}_out.clean.id" \
    | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t'     "${opt_biofile}" - \
    | awk -F $'\t' -v hvm=${cfg_hvm} '{
      OFS="\t"
      hvmtag="OK"
      if ( hvm != 1 ) hvmtag="__NA__"
      print( $0, hvmtag )
    }'
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.2.bio"
mv "${tmpprefix}.2.bio" "${opt_biofile}"

# update biography file with sample relationships
printf "> updating biography file with relatedness information..\n"
{
  extract_related_lists_from_grm_file "${tmpprefix}_sq.genome.gz" \
    | join -t $'\t'     "${opt_biofile}" -
  extract_related_lists_from_grm_file "${tmpprefix}_sq.genome.gz" \
    | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' '{ OFS="\t"; print( $0, "__NA__" ) }'
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.3.bio"
mv "${tmpprefix}.3.bio" "${opt_biofile}"

rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out"*

rm "${tmpprefix}"*

