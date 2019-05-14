#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

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

declare plinkflag=''
# find potential mixups?
if [ ${cfg_hvm} -eq 1 ] ; then
  printf "> computing individual heterozygosity and missing rates..\n"
  ${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}" \
               --set-hh-missing \
               --het \
               --out "${tmpprefix}_sq" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  ${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}" \
               --set-hh-missing \
               --missing \
               --out "${tmpprefix}_sq" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  printf "> removing mixups and low-coverage individuals..\n"
  join --header -t $'\t' \
    <( attach_uids -h "${tmpprefix}_sq.imiss" \
      | tail -n +2               | sort -t $'\t' -u -k 1,1 ) \
    <( attach_uids -h "${tmpprefix}_sq.het" \
      | tail -n +2 | cut -f 1,4- | sort -t $'\t' -u -k 1,1 ) \
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
  plinkflag="--keep ${tmpprefix}_out.clean.id"
else
  # write a dummy clean.id file including everyone
  cut -f 1,2 "${opt_hqprefix}.fam" | sort -u > "${tmpprefix}_out.clean.id"
fi
# extract high coverage variants
tmp_varmiss=${cfg_varmiss}
M=$( wc -l "${opt_inprefix}.fam" | cut -d ' ' -f 1 )
if [ $M -lt ${cfg_minindcount} ] ; then tmp_varmiss=0.1 ; fi
printf "> extracting high coverage variants..\n"
echo "  ${plinkexec} --allow-extra-chr --bfile ${opt_inprefix} ${plinkflag}
             --geno ${tmp_varmiss}
             --make-just-bim
             --out ${tmpprefix}_hcv" | printlog 2
${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" ${plinkflag} \
             --geno ${tmp_varmiss} \
             --make-just-bim \
             --out "${tmpprefix}_hcv" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
cut -d ' ' -f 2 "${tmpprefix}_hcv.bim" > "${tmpprefix}_hcv.mrk"
# extract high coverage individuals
tmp_samplemiss=${cfg_samplemiss}
N=$( wc -l "${opt_inprefix}.bim" | cut -d ' ' -f 1 )
if [ $N -lt ${cfg_minvarcount} ] ; then tmp_samplemiss=0.1 ; fi
printf "> extracting high coverage individuals..\n"
echo "  ${plinkexec} --allow-extra-chr --bfile ${opt_inprefix} ${plinkflag}
             --extract ${tmpprefix}_hcv.mrk
             --mind ${tmp_samplemiss}
             --make-just-fam
             --out ${tmpprefix}_hci" | printlog 2
${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" ${plinkflag} \
             --extract "${tmpprefix}_hcv.mrk" \
             --mind ${tmp_samplemiss} \
             --make-just-fam \
             --out "${tmpprefix}_hci" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# identify related individuals
printf "> identifying related individuals..\n"
${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}" \
             --keep "${tmpprefix}_hci.fam" \
             --set-hh-missing \
             --genome gz \
             --out "${tmpprefix}_sq" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
${plinkexec} --allow-extra-chr --bfile "${opt_hqprefix}" \
             --keep "${tmpprefix}_hci.fam" \
             --set-hh-missing \
             --cluster \
             --read-genome "${tmpprefix}_sq.genome.gz" \
             --rel-cutoff ${cfg_pihatrel} \
             --out "${tmpprefix}_sq" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
# rename list of clean unrelated individuals for later use
mv "${tmpprefix}_sq.rel.id" "${opt_outprefixbase}.ids"
# erase eventual dysfunctional family information
cut -f 1 "${opt_hqprefix}.fam" | sort | uniq -c | tabulate \
  | awk -F $'\t' '{ OFS="\t"; if ( $1 > 2 ) print( $2 ); }' \
  | join -t $'\t' - <( sort -k 1,1 "${opt_hqprefix}.fam" ) \
  | awk '{
      OFS="\t"
      dad[$1,$2] = $1 SUBSEP $3
      mom[$1,$2] = $1 SUBSEP $4
      sex[$1,$2] = $5
    } END{
      for ( sid in dad ) {
        dysfamvec = sex[dad[sid]] == 0 || sex[mom[sid]] == 0
        dysfamvec = dysfamvec || sex[dad[sid]] == sex[mom[sid]]
        if ( dysfamvec ) {
          split( sid, vid, SUBSEP )
          print( vid[1], vid[2], 0, 0 )
        }
      }
    }' \
  > "${tmpprefix}_dysfam.tri"
# remove mixups and update sex and parents in input set
${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" ${plinkflag} \
             --update-parents "${tmpprefix}_dysfam.tri" \
             --update-sex "${opt_hqprefix}.fam" 3 \
             --make-bed \
             --out "${tmpprefix}_out" \
             2> >( tee "${tmpprefix}.err" ) | printlog 3
if [ $? -ne 0 ] ; then
  cat "${tmpprefix}.err"
fi
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"
unset plinkflag

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
{
  {
    printf "%s\t" ${cfg_uid}
    head -n 1 "${opt_hqprefix}.sexcheck" | tabulate | cut -f 3-
  } | join -t $'\t'     "${opt_biofile}" - \
    | tee "${tmpprefix}.0.bio"
  # count number of fields in the merged file
  TNF=$( cat "${tmpprefix}.0.bio" | wc -w )
  # merge information for existing individuals
  attach_uids "${opt_hqprefix}.sexcheck" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t'     "${opt_biofile}" - \
  # add non-existing individuals and pad the extra fields with NAs
  attach_uids "${opt_hqprefix}.sexcheck" \
    | tail -n +2 | cut -f 1,4- | sort -u -k 1,1 \
    | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.1.bio"
mv "${tmpprefix}.1.bio" "${opt_biofile}"

# update biography file with heterozygosity information
{
  {
    printf "%s\t" ${cfg_uid}
    head -n 1 "${tmpprefix}_sq.het" | tabulate | cut -f 3-
  } | join -t $'\t'     "${opt_biofile}" - \
    | tee "${tmpprefix}.2.bio"
  # count number of fields in the merged file
  TNF=$( cat "${tmpprefix}.2.bio" | wc -w )
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
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.3.bio"
mv "${tmpprefix}.3.bio" "${opt_biofile}"

# update biography file with potential mixup information
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
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.4.bio"
mv "${tmpprefix}.4.bio" "${opt_biofile}"

# update biography file with sample relationship
{
  extract_related_lists_from_grm_file "${tmpprefix}_sq.genome.gz" \
    | join -t $'\t'     "${opt_biofile}" -
  extract_related_lists_from_grm_file "${tmpprefix}_sq.genome.gz" \
    | join -t $'\t' -v1 "${opt_biofile}" - \
    | awk -F $'\t' '{ OFS="\t"; print( $0, "__NA__" ) }'
} | sort -t $'\t' -u -k 1,1 > "${tmpprefix}.5.bio"
mv "${tmpprefix}.5.bio" "${opt_biofile}"

rename "${tmpprefix}_out" "${opt_outprefix}" "${tmpprefix}_out".*

rm "${tmpprefix}"*

