#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

bcftoolsexec=${BASEDIR}/lib/3rd/bcftools

declare -r tmpprefix=${opt_outprefix}_tmp
declare -r cfg_chromosomes=$( cfgvar_get chromosomes )

#-------------------------------------------------------------------------------

# input: final QC plink file set
# output: chromosome-wise eagle VCF input files

echo -e "==== Convert ====\n" | printlog 0

printf "\
  * Purge individual and family IDs
  * Convert plink file set into VCF format
\n" | printlog 0

declare -ra cfg_chromosomes_arr=( $cfg_chromosomes )
tmp_chromosomes=$( join \
  <( printf "%s\n" ${cfg_chromosomes_arr[@]} | sort -u ) \
  <( cut -f 1 ${opt_inprefix}.bim | sort -u ) \
)

declare cskip=1
for chr in ${tmp_chromosomes} ; do
  [ -s ${opt_outprefix}_chr${chr}.bcf ] || cskip=0
done
if [ ${cskip} -eq 1 ] ; then
  printf "> all bcf files present. skipping conversion..\n\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "> temporary files '%s*' found. please remove them before re-run.\n\n" "${tmpprefix}" >&2
  exit 1
fi

# purge fam file ids of any unwanted characters
touch ${tmpprefix}_updateids.txt
touch ${tmpprefix}_updateparents.txt
awk -f ${BASEDIR}/lib/awk/idclean.awk \
  -v idmapfile=${tmpprefix}.idmap \
  -v idsupfile=${tmpprefix}_updateids.txt \
  -v parupfile=${tmpprefix}_updateparents.txt \
  --source 'BEGIN{
    fn=0
  } {
    if ( FNR==1 ) fn++
    fid=idclean($1)
    iid=idclean($2)
    gid=idclean($3)
    jid=idclean($4)
    uid=idclean($1"_"$2)
    gsub( "_+$", "", fid )
    gsub( "^_+", "", fid )
    gsub( "_+$", "", iid )
    gsub( "^_+", "", iid )
    gsub( "_+$", "", gid )
    gsub( "^_+", "", gid )
    gsub( "_+$", "", jid )
    gsub( "^_+", "", jid )
    gsub( "_+$", "", uid )
    gsub( "^_+", "", uid )
    if ( fn==1 ) print( $1, $2, fid, iid ) > idmapfile
    if ( fn==2 ) print( uid, uid, gid, jid ) > idsupfile
    if ( fn==3 ) print( fid, iid, gid, jid ) > parupfile
  }' \
    ${opt_inprefix}.fam \
    ${opt_sqprefix}_updateids.txt \
    ${opt_sqprefix}_updateparents.txt
Nold=$( sort -u -k 1,2 ${tmpprefix}.idmap | wc -l )
Nnew=$( sort -u -k 3,4 ${tmpprefix}.idmap | wc -l )
if [ $Nold -eq $Nnew ] ; then
  mv ${tmpprefix}_updateids.txt ${opt_outprefix}_updateids.txt
  mv ${tmpprefix}_updateparents.txt ${opt_outprefix}_updateparents.txt
  echo -e "  ${plinkexec##*/} --allow-extra-chr --bfile ${opt_inprefix}
               --update-ids ${tmpprefix}.idmap
               --make-bed --out ${tmpprefix}_idfix\n" | printlog 2
  ${plinkexec} --allow-extra-chr --bfile "${opt_inprefix}" \
               --update-ids "${tmpprefix}.idmap" \
               --make-bed --out "${tmpprefix}_idfix" \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
else
  printf "====> warning: could not purge IDs due to ensuing conflicts.\n\n"
  awk '{ uid=$1"_"$2; print( uid, uid, $3, $4 ); }' \
    ${opt_sqprefix}_updateids.txt > ${opt_outprefix}_updateids.txt
  cp ${opt_sqprefix}_updateparents.txt ${opt_outprefix}_updateparents.txt
  cp ${opt_inprefix}.bed ${tmpprefix}_idfix.bed
  cp ${opt_inprefix}.bim ${tmpprefix}_idfix.bim
  cp ${opt_inprefix}.fam ${tmpprefix}_idfix.fam
fi
awk '{ print( $1"_"$2, $1"_"$2, $5 ); }' \
  ${tmpprefix}_idfix.fam > ${opt_outprefix}_updatesex.txt
# convert input files to vcf formats
for chr in ${tmp_chromosomes} ; do
  if [ -e ${opt_outprefix}_chr${chr}.bcf ] ; then
    printf "> bcf file '%s' found. skipping conversion..\n\n" ${opt_outprefix}_chr${chr}.bcf
    continue
  fi
  ${plinkexec} --allow-extra-chr --bfile ${tmpprefix}_idfix \
               --chr ${chr} \
               --recode vcf bgz \
               --out ${tmpprefix}_chr${chr} \
               2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  ${bcftoolsexec} convert -Ob --threads 3 \
                  ${tmpprefix}_chr${chr}.vcf.gz \
                  > ${tmpprefix}_chr${chr}_out.bcf.gz
  ${bcftoolsexec} annotate --rename-chrs <( echo -e '23 X\n25 X' ) -Ob --threads 3 \
                  ${tmpprefix}_chr${chr}_out.bcf.gz \
                  > ${tmpprefix}_chr${chr}_outx.bcf.gz
  ${bcftoolsexec} index ${tmpprefix}_chr${chr}_outx.bcf.gz --threads 3
  mv ${tmpprefix}_chr${chr}_outx.bcf.gz     ${opt_outprefix}_chr${chr}.bcf.gz
  mv ${tmpprefix}_chr${chr}_outx.bcf.gz.csi ${opt_outprefix}_chr${chr}.bcf.gz.csi
done

rm -r ${tmpprefix}*

