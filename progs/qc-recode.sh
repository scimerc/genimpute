#!/usr/bin/env bash

# exit on error
trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -ra batchfiles=( ${opt_inputfiles} )

#-------------------------------------------------------------------------------

# for every batch
  # extract var accoding to var whitelist (if enabled)
  # extract samples according to sample whitelist (if enabled)
  # convert colocalized variant set to human readable format (tped) (for later qc)
  # convert to binary plink (for easy use downstream)


for i in ${!batchfiles[@]} ; do
  declare batchcode=$( get_unique_filename_from_path ${batchfiles[$i]} )
  declare boutprefix=${opt_outprefix}_${batchcode}
  declare tmpprefix=${boutprefix}_tmp
  declare debuglogfn=${tmpprefix}_debug.log
  # check for hash collisions
  if [ -f "${boutprefix}.bed" -a -f "${boutprefix}.bim" -a -f "${boutprefix}.fam" ]; then
    printf "'%s' already exists - skipping recode step..\n" "${boutprefix}.bed"
    printf "try increasing 'numchars' in the hash function if you think this should not happen.\n"
    continue
  fi
  if ls ${tmpprefix}* > /dev/null 2>&1; then
    printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
    exit 1
  fi
  echo "converting batch ${batchfiles[$i]}.."
  # define input specific plink settings
  declare flagkeep=''
  declare flagformat='--bfile'
  declare plinkinputfn=${batchfiles[$i]}
  case "$( get_genotype_file_format "${batchfiles[$i]}" )" in
    "bed" )
      flagformat='--bfile'
      plinkinputfn=${batchfiles[$i]%.bed}
      ;;
    "bcf" )
      flagformat='--bcf'
      ;;
    "vcf" )
      flagformat='--vcf'
      ;;
    * )
      printf "error: fileformat not handled\n" >&2
      exit 1
      ;;
  esac
  # use whitelist if existing
  if [ ! -z "${opt_samplewhitelist}" ] ; then
    flagkeep="--keep ${opt_samplewhitelist}"
  fi
  declare bedflag="${flagkeep}"
  declare pedflag="${flagkeep}"
  if [ -s "${opt_varwhitelist}" ] ; then
    bedflag="${bedflag} --extract range ${opt_varwhitelist}"
  fi
  # convert to plink binary format
  plink $flagformat ${plinkinputfn} \
    --merge-x no-fail \
    --make-bed \
    --out ${tmpprefix}_mx \
    2>&1 >> ${debuglogfn} \
    | tee -a ${debuglogfn}
  if [ -z "${bedflag}" ] ; then
    plink $flagformat ${tmpprefix}_mx ${bedflag} \
      --make-bed \
      --out ${tmpprefix}_out \
      2>&1 >> ${debuglogfn} \
      | tee -a ${debuglogfn}
  else
    mv ${tmpprefix}_mx.bed ${tmpprefix}_out.bed
    mv ${tmpprefix}_mx.bim ${tmpprefix}_out.bim
    mv ${tmpprefix}_mx.fam ${tmpprefix}_out.fam
  fi
  # rename variants to universal code chr:bp_a1_a2
  make_variant_names_universal_in_bim_file ${tmpprefix}_out.bim
  # extract colocalized variant positions from gp file and make a tped file set from them
  awk '{ print( $2, $1, 0, $4 ); }' ${tmpprefix}_out.bim \
    | sort -k 1,1 > ${tmpprefix}.gp
  get_variant_info_for_dup_chr_cm_bp_aa_mm ${tmpprefix}.gp \
    | awk '{
        OFS="\t"
        pos0=$4>0?($4-1):0
        print( $2, pos0, $4, $1 )
      }' \
    | sort -u \
    > ${tmpprefix}.coloc.rng
  if [ -s "${tmpprefix}.coloc.rng" ] ; then
    pedflag="${pedflag} --extract range ${tmpprefix}.coloc.rng"
    plink --bfile ${tmpprefix}_out ${pedflag} \
      --recode transpose \
      --out ${tmpprefix}_out \
      2>&1 >> ${debuglogfn} \
      | tee -a ${debuglogfn}
  else
    echo "no colocalized variants found."
    echo "skipping batch '${batchfiles[$i]}' tped recoding.."
    touch ${tmpprefix}_out.tped ${tmpprefix}_out.tfam
  fi
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
  # save fam files for later batch effect dectection
  mv ${tmpprefix}_out.bed ${boutprefix}.bed
  mv ${tmpprefix}_out.bim ${boutprefix}.bim
  mv ${tmpprefix}_out.fam ${boutprefix}.fam
  mv ${tmpprefix}_out.tped ${boutprefix}.tped
  mv ${tmpprefix}_out.tfam ${boutprefix}.tfam
  rm ${tmpprefix}*
  unset flagkeep
  unset flagformat
  unset plinkinputfn
  unset batchcode
  unset boutprefix
  unset tmpprefix
  unset debuglogfn
done

rm ${tmpprefix}*

