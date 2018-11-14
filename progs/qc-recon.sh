#!/usr/bin/env bash

trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -r  tmpprefix=${opt_outprefix}_tmp
declare -r  debuglogfn=${tmpprefix}_debug.log
declare -ra batchfiles=( ${opt_inputfiles} )

declare -r cfg_genomebuild="$( cfgvar_get genomebuild )"
declare -r cfg_refallelesfn="$( cfgvar_get refallelesfn )"

#-------------------------------------------------------------------------------

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping aligment..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
  exit 1
fi


#-------------------------------------------------------------------------------

# for every batch
  # extract var accoding to var whitelist (if enabled)
  # extract samples according to sample whitelist (if enabled)
  # convert colocalized variant set to human readable (for later qc)
  # convert to binary plink (for easy use downstream)

for i in ${!batchfiles[@]} ; do
  echo "converting batch ${batchfiles[$i]}.."
  # define input specific plink settings
  declare flagextract=''
  declare flagformat='--bfile'
  declare plinkinputfn=${batchfiles[$i]}
  declare plinkoutputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_filtered
  declare plinktmpfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_unfiltered
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
    flagextract="--keep ${opt_samplewhitelist}"
  fi
  declare bedflag="${flagextract}"
  declare pedflag="${flagextract}"
  if [ ! -z "${opt_varwhitelist}" ] ; then
    bedflag="${bedflag} --extract range ${opt_varwhitelist}"
  fi
  # check for hash collisions
  if [ -f "${plinkoutputfn}.bed" ]; then
    printf "error: %s already exists - increase 'numchars' in hash function\n" \
      "${plinkoutputfn}.bed" >&2
    exit 1
  fi
  # convert to plink binary format
  plink $flagformat ${plinkinputfn} \
    --merge-x no-fail \
    --make-bed \
    --out ${plinktmpfn} \
    2>&1 >> ${debuglogfn} \
    | tee -a ${debuglogfn}
  if [ -z "${bedflag}" ] ; then
    plink $flagformat ${plinktmpfn} ${bedflag} \
      --make-bed \
      --out ${plinkoutputfn} \
      2>&1 >> ${debuglogfn} \
      | tee -a ${debuglogfn}
  else
    mv ${plinktmpfn}.bed ${plinkoutputfn}.bed
    mv ${plinktmpfn}.bim ${plinkoutputfn}.bim
    mv ${plinktmpfn}.fam ${plinkoutputfn}.fam
  fi
  # rename variants to universal code chr:bp_a1_a2
  make_variant_names_universal_in_bim_file ${plinkoutputfn}.bim
  # extract colocalized variant positions from gp file and make a tped file set from them
  touch ${plinkoutputfn}.tped ${plinkoutputfn}.tfam
  awk '{ print( $2, $1, 0, $4 ); }' ${plinkoutputfn}.bim \
    | sort -k 1,1 > ${plinktmpfn}.gp
  get_variant_info_for_dup_chr_cm_bp_aa ${plinktmpfn}.gp \
    | awk '{
        OFS="\t"
        pos0=$4>0?($4-1):0
        print( $2, pos0, $4, $1 )
      }' \
    | sort -u \
    > ${plinktmpfn}.coloc.rng
  if [ -s "${plinktmpfn}.coloc.rng" ] ; then
    pedflag="${pedflag} --extract range ${plinktmpfn}.coloc.rng"
    plink --bfile ${plinkoutputfn} ${pedflag} \
      --recode transpose \
      --out ${plinkoutputfn} \
      2>&1 >> ${debuglogfn} \
      | tee -a ${debuglogfn}
  else
    echo "no colocalized variants found."
    echo "skipping batch '${batchfiles[$i]}' tped recoding.."
  fi
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${plinkoutputfn}.bim
  sed -i -r 's/[ \t]+/\t/g' ${plinkoutputfn}.fam
  # save fam files for later batch effect dectection
  cp ${plinkoutputfn}.fam \
    ${opt_batchoutprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} ).fam
  unset batchtped
  unset flagextract
  unset flagformat
  unset plinkinputfn
done

