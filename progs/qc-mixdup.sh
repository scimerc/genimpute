#!/usr/bin/env bash

# exit on error
trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log

declare -r cfg_hvm=$( cfgvar_get hvm )
declare -r cfg_pihat=$( cfgvar_get pihat )
declare -r cfg_uid=$( cfgvar_get uid )

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping dup&mix step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
  exit 1
fi

# input: merged plink set and hq plink set
# output: clean plink set (with imputed sex from hq set, no duplicates and no mixups)


# update biography file with sex information
{
  # merge information for existing individuals
  paste_sample_ids ${opt_hqprefix}.sexcheck \
    | join -t $'\t' ${opt_biofile} - \
    | tee ${tmpprefix}.0.bio
  # count number of fields in the merged file
  TNF=$( head -n 1 ${tmpprefix}.0.bio | wc -w )
  # add non-existing individuals and pad the extra fields with NAs
  paste_sample_ids ${opt_hqprefix}.sexcheck \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' -v TNF=${TNF} '{
      OFS="\t"
      printf($0)
      for ( k=NF; k<TNF; k++ ) printf("\t__NA__")
      printf("\n")
    }'
} | sort -u -k 1,1 > ${tmpprefix}.1.bio
cp ${tmpprefix}.1.bio ${opt_biofile}

declare plinkflag=''
# run 'het_VS_miss.Rscript' to find potential mixups?
if [ ${cfg_hvm} -eq 1 ] ; then
  echo "computing individual heterozygosity and missing rates.."
  plink --bfile ${opt_hqprefix} --set-hh-missing --het     --out ${tmpprefix}_sq >> ${debuglogfn}
  plink --bfile ${opt_hqprefix} --set-hh-missing --missing --out ${tmpprefix}_sq >> ${debuglogfn}
  ${BASEDIR}/progs/het_vs_miss.Rscript -m ${tmpprefix}_sq.imiss -h ${tmpprefix}_sq.het \
    -o ${tmpprefix}_out >> ${debuglogfn}
  mv ${tmpprefix}_out_hetVSmiss.pdf ${opt_outprefix}_hetVSmiss.pdf
  mv ${tmpprefix}_out_hetVSmiss.log ${opt_outprefix}_hetVSmiss.log
  [ -s "${tmpprefix}_out.clean.id" ] || {
    printf "error: file '%s' empty or not found.\n" ${tmpprefix}_out.clean.id >&2;
    exit 1;
  }
  # update biography file with potential mixup information
  {
    cut -f 3 ${tmpprefix}_out.clean.id \
      | sort -u \
      | join -t $'\t' -v1 ${opt_biofile} - \
      | awk -F $'\t' '{
        OFS="\t"
        if ( NR == 1 ) print( $0, "MISMIX" )
        else print( $0, "PROBLEM" )
      }'
    cut -f 3 ${tmpprefix}_out.clean.id \
      | sort -u \
      | join -t $'\t'     ${opt_biofile} - \
      | awk -F $'\t' '{
        OFS="\t"
        print( $0, "OK" )
      }'
  } | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.2.bio
  cp ${tmpprefix}.2.bio ${opt_biofile}
  # include non-mixup info later
  plinkflag="--keep ${tmpprefix}_out.clean.id"
fi
# identify duplicate individuals
echo "identifying duplicate individuals.."
plink --bfile ${opt_hqprefix} ${plinkflag} \
      --set-hh-missing \
      --genome gz \
      --out ${tmpprefix}_sq \
      >> ${debuglogfn}
plink --bfile ${opt_hqprefix} ${plinkflag} \
      --set-hh-missing \
      --cluster \
      --read-genome ${tmpprefix}_sq.genome.gz \
      --rel-cutoff ${cfg_pihat} \
      --out ${tmpprefix}_sq \
      >> ${debuglogfn}
# give rel.id file a less confusing name
mv ${tmpprefix}_sq.rel.id ${tmpprefix}_sq.ind
unset plinkflag

extract_related_lists_from_grm_file() {
  local -r infile="$1"
  zcat -f "${infile}" | tabulate \
    | awk -F $'\t' -v uid=${cfg_uid} '{
      OFS="\t"
      maxcnt=11111
      uid0=$1"_"$2
      uid1=$3"_"$4
      if ( NR>1 && $10>=0.1 ) {
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
      printf( "%s\tRELSHIP\n", uid )
      for ( uid in relarr ) print( uid, relarr[uid] )
    }' \
  | sort -t $'\t' -u -k 1,1
}

# update biography file with sample relationship
{
  extract_related_lists_from_grm_file ${tmpprefix}_sq.genome.gz \
    | join -t $'\t' ${opt_biofile} -
  extract_related_lists_from_grm_file ${tmpprefix}_sq.genome.gz \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk '{
      OFS="\t"
      print( $0, "__NA__" )
    }'
} | sort -u -k 1,1 > ${tmpprefix}.3.bio
cp ${tmpprefix}.3.bio ${opt_biofile}

# update biography file with duplicates
{
  awk -F $'\t' '{ print( $1"\t"$2 ); }' ${tmpprefix}_sq.ind \
    | sort -u \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' '{
        OFS="\t"
        if ( NR == 1 ) print( $0, "DUPORG" )
        else {
          if ( $(NF-1) == "PROBLEM" ) print( $0, 0 )
          else print( $0, "DUP" )
        }
      }'
  awk -F $'\t' '{ print( $1"\t"$2 ); }' ${tmpprefix}_sq.ind \
    | sort -u \
    | join -t $'\t' ${opt_biofile} - \
    | awk -F $'\t' '{
        OFS="\t"
        print( $0, "ORG" )
      }'
} | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.4.bio
cp ${tmpprefix}.4.bio ${opt_biofile}

# remove duplicates + update sex in input set
echo "removing duplicate individuals and updating sex.."
plink --bfile ${opt_inprefix} \
      --update-sex ${opt_hqprefix}.fam 3 \
      --keep ${tmpprefix}_sq.ind \
      --make-bed \
      --out ${tmpprefix}_out \
      >> ${debuglogfn}
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam

mv ${tmpprefix}_out.bed ${opt_outprefix}.bed
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.fam ${opt_outprefix}.fam

rm ${tmpprefix}*

