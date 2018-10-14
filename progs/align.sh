#!/usr/bin/env bash

trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r  tmpprefix=${opt_outprefix}_tmp
declare -r  debuglogfn=${tmpprefix}_debug.log
declare -ra batchfiles=( ${opt_inputfiles} )

declare -r genomebuild="b37"

if [ -f "${opt_outprefix}.bed" -a -f "${opt_outprefix}.bim" -a -f "${opt_outprefix}.fam" ] ; then
  printf "skipping aligment..\n"
  exit 0
fi

if ls ${tmpprefix}* > /dev/null 2>&1; then
  printf "error: temporary files exist in '%s'. pls remove\n" "${tmpprefix}" >&2
  exit 1
fi

if ! which bedtools > /dev/null 2>&1; then
  printf "error: bedtools not found\n" >&2
  exit 1
fi


get_genotype_file_format()  {
  
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


get_unique_filename_from_path() {
  local -r path="$1"
  local -r numchars=6
  local -r md5=$(echo "$path" | md5sum)
  printf "%s_%s" "$(basename "$path" )" "${md5:0:$numchars}"
}


# in: plink tped; out: (rs, chr, cm, bp) 
get_variant_info_from_tped() {
  awk '{
    delete catalog
    for ( k = 5; k <= NF; k++ ) {
      if ( $k > 0 ) catalog[$k] = 1
    }
    asorti( catalog )
    # write the variant name first to enable skipping it in the uniq command later on
    printf( "%s %s %s %s", $2, $1, $3, $4 ) # rs,chr,cm,bp
    for ( allele in catalog ) printf( "_%s", catalog[allele] )
    print( "" )
  }' -
}


# in: (rs, chr, cm, bpa); out: same
# get unique chr,cm,bpa entries [sort|uniq]
# sort on rsnumbers corresponding to unique entries
# get complementary set of rsnumbers (rsnumbers of duplicated vars) [join -v1 | cut]
get_variant_info_for_dup_chr_cm_bpa() {
  local -r inputfile="$1"
  sort -k 2 ${inputfile} \
    | uniq -u -f 1 \
    | sort -t ' ' -k 1,1 \
    | join -t ' ' -v1 ${inputfile} -
}


# whitelist of variants
declare -r varwhitelist=${tmpprefix}.extract


# compute whitelist?
if [ $opt_mini -eq 1 ] ; then
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
    declare extmp=${tmpprefix}_extmp
    # whatever format the input file is - make a bim file
    plink ${plinkflag} ${batchfiles[$i]/%.bed/.bim} --make-just-bim --out ${extmp} >> ${debuglogfn}
    perl -p -i -e 's/[ \t]+/\t/g' ${extmp}.bim
    if [ -s "${extmp}.bim" ] ; then
      if [ $i -eq 0 ] ; then
        # make the other type (non-plink) of bed file from the bim file
        awk -F $'\t' '{ OFS="\t"; print( $1, $4 - 1, $4, $2 ); }' ${extmp}.bim > ${extmp}.bed
      else
        bedtools intersect \
          -a ${varwhitelist} \
          -b <( awk -F $'\t' '{ OFS="\t"; print( $1, $4-1, $4, $2 ); }' ${extmp}.bim ) \
        > ${extmp}.bed
      fi
      mv ${extmp}.bed ${varwhitelist}
      rm -f ${extmp}.bim
    fi
    unset plinkflag
    unset extmp
  done
fi




# for every batch
  # extract var accoding to var whitelist (if enabled)
  # extract samples according to sample whitelist (if enabled)
  # convert to human readable (for later qc)
  # convert to binary plink (for easy use downstream)

for i in ${!batchfiles[@]} ; do
  echo "converting batch ${batchfiles[$i]}.."
  # define input specific plink settings
  declare flagextract=''
  declare flagformat='--bfile'
  declare plinkinputfn=${batchfiles[$i]}
  declare plinkoutputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_filtered
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
    flagextract="${flagextract} --keep ${opt_samplewhitelist}"
  fi
  if [ ! -z "${varwhitelist}" ] ; then
    flagextract="${flagextract} --extract range ${varwhitelist}"
  fi
  # check for hash collisions
  if [ -f "${plinkoutputfn}.bed" ]; then
    printf "error: %s already exists - increase 'numchars' in hash function\n" \
      "${plinkoutputfn}.bed" >&2
    exit 1
  fi
  {
    plink $flagformat ${plinkinputfn} ${flagextract} --recode transpose --out ${plinkoutputfn}
    plink $flagformat ${plinkinputfn} ${flagextract} --make-bed         --out ${plinkoutputfn}
  } >> ${debuglogfn}
  # tab-separate all human-readable plink files
  perl -p -i -e 's/[ \t]+/\t/g' ${plinkoutputfn}.bim
  perl -p -i -e 's/[ \t]+/\t/g' ${plinkoutputfn}.fam
  # save fam files for batch effect dectection
  cp ${plinkoutputfn}.fam \
    ${opt_batchoutprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} ).fam
  unset batchtped
  unset flagextract
  unset flagformat
  unset plinkinputfn
done


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
  declare plinkinputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_filtered
  declare plinkoutputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_purged
  declare tmpvarctrl=${tmpprefix}_varctrl
  declare tmpvardups=${tmpprefix}_vardups
  declare batchblacklist=${plinkinputfn}.blacklist
  declare batchfliplist=${plinkinputfn}.fliplist
  # extract marker information corresponding to unique position-specific genotype series:
  # identical position-specific genotype series will thereby be ignored as harmless here.
  sort -t ' ' -u -k 1,1 -k 4 ${plinkinputfn}.tped \
    | get_variant_info_from_tped \
    | sort -t ' ' -k 1,1 \
    > ${plinkinputfn}.gp
  get_variant_info_for_dup_chr_cm_bpa ${plinkinputfn}.gp \
    | cut -d ' ' -f 1 \
    | uniq \
    > ${batchblacklist} # format: rs
  #TODO: check actual role of centimorgans here
  echo -n "$( wc -l ${batchblacklist} | cut -d ' ' -f 1 ) "
  echo "non-coherent duplicate variants marked for deletion."
  # extract marker information corresponding to non-blacklisted genotype series
  sort -t ' ' -k 2,2 ${plinkinputfn}.tped \
    | get_variant_info_from_tped \
    | join -t ' ' -v1 - ${batchblacklist} \
    > ${plinkinputfn}.gpz
  # get variant with highest cm from each set of duplicates [sort -k 3,3 | sort -u]
  # get their rsnumbers [cut | sort -u]
  get_variant_info_for_dup_chr_cm_bpa ${plinkinputfn}.gpz \
    | sort -t ' ' -k 3,3r \
    | sort -t ' ' -u -k 2,2 -k 4,4 \
    | cut -d ' ' -f 1 \
    | sort -u \
    > ${tmpvardups}  
  declare wl_size="$( wc -l ${tmpvardups} | cut -d ' ' -f 1 )"
  # add rsnumbers of non-retained duplicated variants to blacklist [join -v1]
  # (count them in the process)
  declare bl_size_old="$( wc -l ${batchblacklist} | cut -d " " -f 1 )"
  get_variant_info_for_dup_chr_cm_bpa ${plinkinputfn}.gpz \
    | join -t ' ' -v1 - ${tmpvardups} \
    >> ${batchblacklist} 
  declare bl_size_new="$( wc -l ${batchblacklist} | cut -d " " -f 1 )"
  rm ${tmpvardups}
  printf "%s unique coherent duplicate variants retained [%s marked for deletion].\n" \
    "${wl_size}" $(( bl_size_new - bl_size_old ))
  # use ref alleles specified or if we have more than one batch
  if [ ${#batchfiles[@]} -gt 1 -o ! -z ${opt_refallelesfn} ] ; then
    echo "matching variants to reference.."
    [ -s ${opt_refallelesfn} ] \
      || { printf "error: file '%s' empty or not found\n" ${opt_refallelesfn} >&2; exit 1; }
    # get chr:bp strings from bim file and join with the corresponding field of refallelesfn 
    awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${plinkinputfn}.bim \
      | sort -t $'\t' -k 7,7 \
      | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' ${opt_refallelesfn} - \
      | awk -F $'\t' \
        -f ${BASEDIR}/lib/awk/nucleocode.awk \
        -f ${BASEDIR}/lib/awk/genotype.awk \
        -f ${BASEDIR}/lib/awk/gflip.awk \
        -f ${BASEDIR}/lib/awk/gmatch.awk \
        -v batchblacklist=${batchblacklist} \
        -v batchfliplist=${batchfliplist} \
        --source 'BEGIN{
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
  if [ -s "${batchblacklist}" ] ; then
    plinkflags="${plinkflags} --exclude ${batchblacklist}"
  fi
  if [ -s "${batchfliplist}" ] ; then
    plinkflags="${plinkflags} --flip ${batchfliplist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $plinkinputfn $plinkoutputfn"
  plink --bfile ${plinkinputfn} ${plinkflags} --make-bed --out ${plinkoutputfn} >> ${debuglogfn}
  # tab-separate all human-readable plink files
  perl -p -i -e 's/[ \t]+/\t/g' ${plinkoutputfn}.bim
  perl -p -i -e 's/[ \t]+/\t/g' ${plinkoutputfn}.fam
  # rename variants to universal code chr:bp_a1_a2
  awk -F $'\t' '{
    OFS="\t";
    a[1] = $5; a[2] = $6; asort(a);
    $2 = $1":"$4"_"a[1]"_"a[2];
    print;
  }' ${plinkoutputfn}.bim > ${tmpvarctrl}
  mv ${tmpvarctrl} ${plinkoutputfn}.bim
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
    declare plinkinputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_purged
    declare plinkoutputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )_ready2merge
    plinkflag=''
    if [ -s "${varblacklist}" ] ; then
      plinkflag="--exclude ${varblacklist}" 
    fi
    plink --bfile ${plinkinputfn} ${plinkflag} --make-bed --out ${plinkoutputfn} >> ${debuglogfn}
    echo "${plinkoutputfn}" >> ${batchlist}
  done 
  # NOTE: if only one batch is present plink throws a warning
  plink --merge-list ${batchlist} --out ${tmpprefix}_out >> ${debuglogfn}
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
  perl -p -i -e 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
  perl -p -i -e 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
  # are we done?
  if [ ! -f "${tmpprefix}.missnp" ] ; then
    break
  fi
  # no! prepare to repeat loop
  mv ${tmpprefix}.missnp ${varblacklist}
  echo "repeating merging attempt.."
done

# pre-process sex chromosomes variants
plinkflag=''
parcount=$( awk '$1 == 25' ${tmpprefix}_out.bim | wc -l )
if [ $parcount -eq 0 ] ; then
  plinkflag="--split-x ${genomebuild} no-fail" 
fi
plink --bfile ${tmpprefix}_out ${plinkflag} --make-bed --out ${tmpprefix}_outsx >> ${debuglogfn}

mv ${tmpprefix}_outsx.bed ${opt_outprefix}.bed
mv ${tmpprefix}_outsx.bim ${opt_outprefix}.bim
mv ${tmpprefix}_outsx.fam ${opt_outprefix}.fam

rm ${tmpprefix}*

