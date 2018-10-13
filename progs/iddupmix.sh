#!/usr/bin/env bash

# exit on error
trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r debuglogfn=${tmpprefix}_debug.log
declare -r cfg_hvm=1
declare -r cfg_pihat=0.9
declare -r cfg_uid=000UID

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping dup&mix step..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
  exit 1
fi

# input: merged plink set and hq plink set
# output: clean plink set (with imputed sex from hq set, no duplicates and no mixups)


declare plinkflag=''
# run 'het_VS_miss.Rscript' to find potential mixups?
if [ ${cfg_hvm} -eq 1 ] ; then
  echo "computing individual heterozygosity and missing rates.."
  plink --bfile ${opt_hqprefix} --het     --out ${tmpprefix}_sq >> ${debuglogfn}
  plink --bfile ${opt_hqprefix} --missing --out ${tmpprefix}_sq >> ${debuglogfn}
  ${BASEDIR}/progs/het_vs_miss.Rscript -m ${tmpprefix}_sq.imiss -h ${tmpprefix}_sq.het \
    -o ${tmpprefix} >> ${debuglogfn}
  [ -s "${tmpprefix}.clean.id" ] || { printf "error: in '%s'\n" ${tmpprefix}.clean.id; exit 1; }
  # update biography file with potential mixup information
  {
    cut -f 3 ${tmpprefix}.clean.id | sort -u | join -t $'\t' -v1 ${opt_biofile} - | awk -F $'\t' '{
      OFS="\t"
      if ( NR == 1 ) print( $0, "het_VS_miss" )
      else print( $0, 0 )
    }'
    cut -f 3 ${tmpprefix}.clean.id | sort -u | join -t $'\t'     ${opt_biofile} - | awk -F $'\t' '{
      OFS="\t"
      print( $0, 1 )
    }'
  } | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.tmpbio
  mv ${tmpprefix}.tmpbio ${opt_biofile}
  # include non-mixup info later 
  plinkflag="--keep ${tmpprefix}.clean.id"
fi
echo "identifying duplicate individuals.."
plink --bfile ${opt_hqprefix} ${plinkflag} \
      --genome gz \
      --out ${tmpprefix}_sq \
      >> ${debuglogfn}
plink --bfile ${opt_hqprefix} ${plinkflag} \
      --cluster \
      --read-genome ${tmpprefix}_sq.genome.gz \
      --rel-cutoff ${pihat} \
      --out ${tmpprefix}_sq \
      >> ${debuglogfn}
# update biography file with sample relationship
zcat ${tmpprefix}_sq.genome.gz \
  | sed -r 's/[ \t]+/\t/g; s/^[ \t]+//g;' \
  | awk -F $'\t' -v uid=${uid} 'BEGIN{
      OFS="\t"
      printf( "%s\tRELSHIP\n", uid )
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
    }' \
  | sort -t $'\t' -u -k 1,1 \
  | join -t $'\t' -a1 -e '-' ${opt_biofile} - \
  > ${tmpprefix}.tmpbio
mv ${tmpprefix}.tmpbio ${opt_biofile}
# update biography file with duplicates
{
  awk -F $'\t' '{ print( $1"\t"$2 ); }' ${tmpprefix}_sq.rel.id \
    | sort -u \
    | join -t $'\t' -v1 ${opt_biofile} - \
    | awk -F $'\t' '{
        OFS="\t"
        if ( NR == 1 ) print( $0, "duplicate" )
        else print( $0, 1 )
      }'
  awk -F $'\t' '{ print( $1"\t"$2 ); }' ${tmpprefix}_sq.rel.id \
    | sort -u \
    | join -t $'\t' ${opt_biofile} - \
    | awk -F $'\t' '{
        OFS="\t"
        print( $0, 1 )
      }'
} | sort -t $'\t' -u -k 1,1 > ${tmpprefix}.tmpbio
mv ${tmpprefix}.tmpbio ${opt_biofile}
# remove duplicates + update sex in input set
plink --bfile ${opt_inprefix} \
      --update-sex ${opt_hqprefix}.fam 3 \
      --keep ${tmpprefix}_sq.rel.id \
      --make-bed \
      --out ${tmpprefix}_out \
      >> ${debuglogfn}
perl -p -i -e 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
perl -p -i -e 's/[ \t]+/\t/g' ${tmpprefix}_out.fam

mv ${tmpprefix}_out.fam ${opt_outprefix}.fam
mv ${tmpprefix}_out.bim ${opt_outprefix}.bim
mv ${tmpprefix}_out.bed ${opt_outprefix}.bed

rm ${tmpprefix}*

