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

get_genotype_file_format() {
  local -r inputfile=$1
  #TODO implement plink2 native .gen format
  local -r bedhex='0x6c 0x1b'
  local -r bcfhex='0x42 0x43'
  local -r vcfhex='0x23 0x23'

  # get first 2 bytes of inputfile
  if [ -z "${inputfile}" -o ! -f "${inputfile}" ]; then
    return 1
  fi

  local -r filehead=$( zcat -f ${inputfile} | hexdump -n2 -e '2/1 "0x%02x "' )
  case ${filehead} in
    "${bedhex}" )
      printf 'bed' ;;
    "${bcfhex}" )
      printf 'bcf' ;;
    "${vcfhex}" )
      printf 'vcf' ;;
    * )
      printf "error: unknown file format '%s'\n" "${filehead}" >&2
      return 1 ;;
  esac
  return 0
}

#-------------------------------------------------------------------------------

get_unique_filename_from_path() {
  local -r path="$1"
  local -r numchars=6
  local -r md5=$(echo "$path" | md5sum)
  printf "%s_%s" "$(basename "$path" )" "${md5:0:$numchars}"
}

#-------------------------------------------------------------------------------

# stdin: plink tped; stdout: (rs, chr, cm, bp)
get_variant_info_from_tped() {
  awk '{
    delete catalog
    for ( k = 5; k <= NF; k++ ) {
      if ( $k > 0 ) catalog[$k] = 1
    }
    asorti( catalog )
    # write the variant name first to enable skipping it in the uniq command later on
    printf( "%s %s %s %s ", $2, $1, $3, $4 ) # rs,chr,cm,bp
    for ( h = 1; h <= length( catalog ); h++ ) {
      if ( h > 1 ) printf( "," )
      printf( "%s", catalog[h] )
    }
    print( "" )
  }' -
}

#-------------------------------------------------------------------------------

# in: (rs, chr, cm, bp[, a]) sorted on rs; stdout: same
# get unique chr,cm,bp[,a] entries [sort|uniq]
# sort on rsnumbers corresponding to unique entries
# get complementary set of entries (duplicated vars) [join -v1]
get_variant_info_for_dup_chr_cm_bp_aa() {
  local -r inputfile="$1"
  sort -k 2 ${inputfile} \
    | uniq -u -f 1 \
    | sort -t ' ' -k 1,1 \
    | join -t ' ' -v1 ${inputfile} -
}

#-------------------------------------------------------------------------------

# in: tab-separated plink bim; out: same
make_variant_names_universal_in_bim_file() {
  local -r fn=$1
  awk -F $'\t' '{
    OFS="\t";
    a[1] = $5; a[2] = $6; asort(a);
    $2 = $1":"$4"_"a[1]"_"a[2];
    print;
  }' ${fn} > ${fn}.tmp
  mv ${fn}.tmp ${fn}
}

#-------------------------------------------------------------------------------

# in: tab-separated plink bim; stdout: chr:bp, rs
bimtogprs() {
  local -r inputfile="$1"
  awk -F $'\t' '{
    OFS="\t"
    print( $1":"$4, $2 )
  }' ${inputfile} \
  | sort -u -k 1,1
}

#-------------------------------------------------------------------------------

# whitelist of variants
declare -r varwhitelist=${tmpprefix}.extract

# compute whitelist?
if [ $opt_minivarset -eq 1 ] ; then
  for i in ${!batchfiles[@]} ; do
    declare plinkflag=""
    case "$( get_genotype_file_format "${batchfiles[$i]}" )" in
      "bed" ) 
        plinkflag="--bim"
        ;;
      "vcf" )
        plinkflag="--vcf"
        ;;
      "bcf" )
        plinkflag="--bcf"
        ;;
      * ) 
        printf "error: fileformat not handled\n" >&2
        exit 1
        ;;
    esac
    # whatever format the input file is - make a bim file
    plink ${plinkflag} ${batchfiles[$i]/%.bed/.bim} \
          --make-just-bim \
          --out ${tmpprefix}_ex \
          2>&1 >> ${debuglogfn} \
          | tee -a ${debuglogfn}
    sed -i -r 's/[ \t]+/\t/g; s/^chr//g;' ${tmpprefix}_ex.bim
    sed -i -r 's/^XY/X/g; s/^X/23/g; s/^Y/24/g; s/^25/23/g;' ${tmpprefix}_ex.bim
    if [ -s "${tmpprefix}_ex.bim" ] ; then
      if [ $i -eq 0 ] ; then
        bimtogprs ${tmpprefix}_ex.bim \
          | sort -t $'\t' -u -k 1,1 \
          > ${tmpprefix}_ex.0.gprs
      else
        join -t $'\t' ${tmpprefix}_ex.1.gprs <( bimtogprs ${tmpprefix}_ex.bim ) \
          | cut -f 1,2 \
          | sort -t $'\t' -u -k 1,1 \
          > ${tmpprefix}_ex.0.gprs
      fi
      mv ${tmpprefix}_ex.0.gprs ${tmpprefix}_ex.1.gprs
      rm -f ${tmpprefix}_ex.bim
    fi
    awk -F $'\t' '{
      OFS="\t"
      split( $1, gpvec, ":" )
      print( gpvec[1], gpvec[2] - 1, gpvec[2], $2 )
    }' ${tmpprefix}_ex.1.gprs > ${varwhitelist}
    unset plinkflag
  done
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
  if [ ! -z "${varwhitelist}" ] ; then
    bedflag="${bedflag} --extract range ${varwhitelist}"
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


#-------------------------------------------------------------------------------


# for every batch
  # create a batch-chr-specific blacklist with non-coherent variants
  # append non compatible ref variants to blacklist
  # create temp fliplist
  # purge batchfile of blacklist
  # align fliplist in batchfile
# keep blacklists and fliplists just in case

for i in ${!batchfiles[@]} ; do
  echo "purging and aligning batch ${batchfiles[$i]}.."
  # define input specific plink settings
  declare tmpvarctrl=${tmpprefix}_varctrl
  declare plinkinputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_filtered
  declare plinkoutputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_aligned
  declare plinktmpfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_purged
  declare batchblacklist=${plinkinputfn}.blacklist
  declare batchfliplist=${plinkinputfn}.fliplist
  # extract marker information corresponding to unique position-specific genotype series:
  # identical position-specific genotype series will thereby be ignored as harmless here.
  sort -t ' ' -u -k 1,1 -k 4 ${plinkinputfn}.tped \
    | get_variant_info_from_tped \
    | sort -t ' ' -k 1,1 \
    > ${plinkinputfn}.gp
  # get variant info from each set of duplicates
  # get their rsnumbers [cut | sort -u]
  #TODO: check role of centimorgans
  get_variant_info_for_dup_chr_cm_bp_aa ${plinkinputfn}.gp \
    | cut -d ' ' -f 1 \
    | sort -u \
    > ${batchblacklist} # format: rs
  # extract marker information corresponding to non-blacklisted genotype series
  sort -k 2,2 ${plinkinputfn}.bim \
    | join -t $'\t' -v2 -2 2 -o '2.1 2.2 2.3 2.4 2.5 2.6' ${batchblacklist} - \
    | awk -F $'\t' -v batchblacklist=${batchblacklist} '{
        OFS="\t"
        if ( $2 in catalog ) {
          catalog[$2]++
          varname = $2"_"catalog[$2]
          if ( $3 > cmcatalog[$2] ) {
            cmcatalog[$2] = $3
            print $2 >>batchblacklist
          } else print varname >>batchblacklist
          $2 = varname
        } else {
          catalog[$2] = 1
          cmcatalog[$2] = $3
        }
        print
      }' ${plinkinputfn}.bim \
    > ${tmpvarctrl}
  mv ${tmpvarctrl} ${plinkinputfn}.bim
  echo -n "$( wc -l ${batchblacklist} | cut -d ' ' -f 1 ) "
  echo "duplicate variants marked for deletion."
  # use ref alleles if specified or if we have more than one batch
  if [ ${#batchfiles[@]} -gt 1 -o ! -z "${cfg_refallelesfn}" ] ; then
    echo "matching variants to reference.."
    [ -s ${cfg_refallelesfn} ] || {
      printf "error: file '%s' empty or not found.\n" "${cfg_refallelesfn}" >&2;
      exit 1;
    }
    # get chr:bp strings from bim file and join with the corresponding field of refallelesfn
    awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${plinkinputfn}.bim \
      | sort -t $'\t' -k 7,7 \
      | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' ${cfg_refallelesfn} - \
      | awk -F $'\t' \
        -f ${BASEDIR}/lib/awk/nucleocode.awk \
        -f ${BASEDIR}/lib/awk/genotype.awk \
        -f ${BASEDIR}/lib/awk/gflip.awk \
        -f ${BASEDIR}/lib/awk/gmatch.awk \
        -v batchblacklist=${batchblacklist} \
        -v batchfliplist=${batchfliplist} \
        --source 'BEGIN{
            OFS="\t"
            total_miss = 0
            total_mism = 0
            total_flip = 0
            printf( "" ) >>batchblacklist
            printf( "" ) >>batchfliplist
          } {
            if ( $5 == "-" || $6 == "-" ) {
              print( $4 ) >>batchblacklist
              total_miss++
            }
            else {
              if ( !gmatch( $2, $3, $5, $6 ) ) {
                print( $4 ) >>batchblacklist
                total_mism++
              }
              else if ( gflip( $2, $3, $5, $6 ) ) {
                print( $4 ) >>batchfliplist
                total_flip++
              }
            }
          } END{
            print( "total missing:  ", total_miss )
            print( "total mismatch: ", total_mism )
            print( "total flipped:  ", total_flip )
            close( batchblacklist )
            close( batchfliplist)
          }'
    # list unique
    sort -u ${batchfliplist} > $tmpvarctrl
    mv $tmpvarctrl ${batchfliplist}
    echo "$( wc -l ${batchfliplist} ) variants to be flipped."
  fi
  # list unique
  sort -u ${batchblacklist} > $tmpvarctrl
  mv $tmpvarctrl ${batchblacklist}
  echo "$( wc -l ${batchblacklist} ) variants to be excluded."

  declare plinkflags=""
  declare plinkflags_eff=""
  if [ -s "${batchblacklist}" ] ; then
    plinkflags_eff="${plinkflags} --exclude ${batchblacklist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $plinkinputfn $plinkoutputfn"
  plink --bfile ${plinkinputfn} ${plinkflags_eff} \
        --make-bed \
        --out ${plinktmpfn} \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  if [ -s "${batchfliplist}" ] ; then
    plinkflags_eff="${plinkflags} --flip ${batchfliplist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $plinkinputfn $plinkoutputfn"
  plink --bfile ${plinktmpfn} ${plinkflags_eff} \
        --make-bed \
        --out ${plinkoutputfn} \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflags
  unset plinkflags_eff
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${plinkoutputfn}.bim
  sed -i -r 's/[ \t]+/\t/g' ${plinkoutputfn}.fam
  # rename variants to universal code chr:bp_a1_a2
  make_variant_names_universal_in_bim_file ${plinkoutputfn}.bim
done


#-------------------------------------------------------------------------------

# init global blacklist
# while true
  # for every batch
    # purge batchfile of global blacklist
  # attempt to merge resulting batchfiles
  # stop on success
  # if blacklist is not empty, tell francesco
  #   (if we run more than 2 times something might be weird)
  # append to global blacklist (derived from errors)


echo "merging batches.."
declare -r varblacklist=${tmpprefix}.exclude
declare -r mismatchlist=${tmpprefix}.mismatch
declare -r batchlist=${tmpprefix}.batchlist
while true ; do
  > ${batchlist}
  for i in ${!batchfiles[@]} ; do
    declare plinkinputfn=${tmpprefix}_$( \
      get_unique_filename_from_path ${batchfiles[$i]}
    )_aligned
    declare plinkoutputfn=${tmpprefix}_$( \
      get_unique_filename_from_path ${batchfiles[$i]}
    )_mergend
    declare plinkflag=''
    if [ -s "${varblacklist}" ] ; then
      plinkflag="--exclude ${varblacklist}" 
    fi
    plink --bfile ${plinkinputfn} ${plinkflag} \
          --make-bed \
          --out ${plinkoutputfn} \
          2>&1 >> ${debuglogfn} \
          | tee -a ${debuglogfn}
    echo "${plinkoutputfn}" >> ${batchlist}
  done
  unset plinkflag
  # NOTE: if only one batch is present plink throws a warning
  plink --merge-list ${batchlist} \
        --out ${tmpprefix}_out \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  # extract plink's warnings about chromosome and position clashes from plink's log and add the
  # corresponding variant names to plink's own missnp file.
  egrep '^Warning: Multiple' ${tmpprefix}_out.log \
    | cut -d ' ' -f 7 \
    | tr -d "'." \
    | sort -u \
    >> ${mismatchlist}
  if [[ -s "${mismatchlist}" ]] ; then
    sort -u ${mismatchlist} >> ${tmpprefix}_out.missnp
  fi
  # are we done?
  if [ ! -f "${tmpprefix}.missnp" ] ; then
    break
  fi
  # no! prepare to repeat loop
  sort -u ${tmpprefix}.missnp > ${varblacklist}
  echo "repeating merging attempt.."
  rm ${tmpprefix}.missnp
done

# pre-process sex chromosomes variants
declare plinkflag=''
parcount=$( awk '$1 == 25' ${tmpprefix}_out.bim | wc -l )
if [ $parcount -eq 0 ] ; then
  plinkflag="--split-x ${cfg_genomebuild} no-fail" 
fi
plink --bfile ${tmpprefix}_out ${plinkflag} \
      --make-bed \
      --out ${tmpprefix}_outsx \
      2>&1 >> ${debuglogfn} \
      | tee -a ${debuglogfn}
unset plinkflag
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_outsx.bim
sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_outsx.fam

mv ${tmpprefix}_outsx.bed ${opt_outprefix}.bed
mv ${tmpprefix}_outsx.bim ${opt_outprefix}.bim
mv ${tmpprefix}_outsx.fam ${opt_outprefix}.fam

rm ${tmpprefix}*

