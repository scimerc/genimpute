#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare  cfg_genomebuild
         cfg_genomebuild="$( cfgvar_get genomebuild )"
readonly cfg_genomebuild

declare -r cfg_samplewhitelist=$( cfgvar_get samplewhitelist )

if [ ! -z "${cfg_refprefix}" ] ; then
  declare -a batchfiles=( ${opt_inputfiles} "${cfg_refprefix}.all.haplotypes.bcf.gz" )
else
  declare -a batchfiles=( ${opt_inputfiles} )
fi
readonly batchfiles

#-------------------------------------------------------------------------------

# input: original genotype file(s)
# output: plink bed and eventual tped file sets

echo -e "==== Recode ====\n" | printlog 0

printf "\
  * For every batch:
    * Extract variants from variant whitelist (if enabled)
    * Extract samples according to sample whitelist (if enabled)
    * Convert colocalized variant set to human readable format (tped) (for later QC)
    * Convert to binary plink (for easy use downstream)
\n" | printlog 0

for i in ${!batchfiles[@]} ; do
  declare batchcode=$( get_unique_filename_from_path "${batchfiles[$i]}" )
  declare b_outprefix="${opt_outprefix}_batch_${batchcode}"
  declare tmpprefix="${b_outprefix}_tmp"
  # check for hash collisions
  if [ -f "${b_outprefix}.bed" -a -f "${b_outprefix}.bim" -a -f "${b_outprefix}.fam" ]; then
    printf "> '%s' already exists. skipping recode step..\n" "${b_outprefix}.bed"
    printf "> increase 'numchars' in the hash function if you think this shouldn't happen.\n\n"
    continue
  fi
  if ls "${tmpprefix}"* > /dev/null 2>&1; then
    printf "> temporary files '%s*' found. please remove them before re-run.\n\n" "${tmpprefix}" >&2
    exit 1
  fi
  printf "> recoding batch '$( basename "${batchfiles[$i]}" )'..\n\n"
  # define input specific plink settings
  declare flagkeep=''
  declare flagformat='--bfile'
  declare plinkinputfn="${batchfiles[$i]}"
  [ -f "${batchfiles[$i]}" ] || { printf "file '%s' not found." "${batchfiles[$i]}"; exit 1; }
  fformat=$( get_genotype_file_format "${batchfiles[$i]}" )
  case "${fformat}" in
    "bed" )
      flagformat='--bfile'
      plinkinputfn="${batchfiles[$i]%.bed}"
      ;;
    "bcf" )
      flagformat='--bcf'
      ;;
    "vcf" )
      flagformat='--vcf'
      ;;
    * )
      printf "> error: unhandled fileformat '%s'.\n\n" ${fformat} >&2
      exit 1
      ;;
  esac
  # use whitelist if existing
  if [ ! -z "${cfg_samplewhitelist}" ] ; then
    flagkeep="--keep ${cfg_samplewhitelist}"
  fi
  declare bedflag=""
  declare pedflag=""
  if [ -s "${opt_varwhitelist}" ] ; then
    bedflag="--extract range ${opt_varwhitelist}"
  fi
  # convert to plink binary format
  echo -e "  ${plinkexec##*/} --allow-extra-chr $flagformat ${plinkinputfn} ${flagkeep}
          --make-bed
          --out ${tmpprefix}_ex\n" | printlog 2
  ${plinkexec} --allow-extra-chr $flagformat "${plinkinputfn}" ${flagkeep} \
    --make-bed \
    --out "${tmpprefix}_ex" \
    2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  # re-write variant info in universal format
  standardize_bim_file "${tmpprefix}_ex.bim"
  parcount=$( awk '$1 == 25' "${tmpprefix}_ex.bim" | wc -l )
  # split X chromosome variants if necessary
  if [ $parcount -eq 0 ] ; then
    echo -e "  ${plinkexec##*/} --allow-extra-chr --bfile ${tmpprefix}_ex
                 --split-x ${cfg_genomebuild} no-fail
                 --make-bed
                 --out ${tmpprefix}_draftex\n" | printlog 2
    ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_ex" \
                 --split-x ${cfg_genomebuild} no-fail \
                 --make-bed \
                 --out "${tmpprefix}_draftex" \
                 2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    rename "${tmpprefix}_draftex" "${tmpprefix}_ex" "${tmpprefix}_draftex".*
  fi
  # extract region whitelist
  echo -e "  ${plinkexec##*/} --allow-extra-chr --bfile ${tmpprefix}_ex ${bedflag}
          --make-bed
          --out ${tmpprefix}_out\n" | printlog 2
  ${plinkexec} --allow-extra-chr --bfile "${tmpprefix}_ex" ${bedflag} \
    --make-bed \
    --out "${tmpprefix}_out" \
    2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  mv "${tmpprefix}_out.log" "${b_outprefix}.1.log"
  # extract colocalized variant positions from gp file and make a tped file set from them
  awk '{ print( $2, $1, 0, $4 ); }' "${tmpprefix}_out.bim" | sort -k 1,1 > "${tmpprefix}.gp"
  get_variant_info_for_dup_chr_cm_bp_aa_mm "${tmpprefix}.gp" \
    | awk '{
        OFS="\t"
        pos1=$4
        pos0=$4>0?($4-1):0
        print( $2, pos0, pos1, $1 )
      }' \
    | sort -u \
    > "${tmpprefix}.coloc.rng"
  if [ -s "${tmpprefix}.coloc.rng" ] ; then
    pedflag="--extract range ${tmpprefix}.coloc.rng"
    ${plinkexec} --allow-extra-chr \
      --bfile "${tmpprefix}_out" ${pedflag} \
      --recode transpose \
      --out "${tmpprefix}_out" \
      2> >( tee "${tmpprefix}.err" ) | printlog 3
    if [ $? -ne 0 ] ; then
      cat "${tmpprefix}.err"
    fi
    mv "${tmpprefix}_out.log" "${b_outprefix}.2.log"
  else
    printf "> no colocalized variants found.\n"
    printf "> skipping batch '$( basename "${batchfiles[$i]}" )' tped recoding..\n\n"
    touch "${tmpprefix}_out.tped" "${tmpprefix}_out.tfam"
  fi
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.bim"
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_out.fam"

  rename "${tmpprefix}_out" "${b_outprefix}" "${tmpprefix}_out".*

  rm "${tmpprefix}"*

  unset flagkeep
  unset flagformat
  unset plinkinputfn
  unset batchcode
  unset b_outprefix
done

