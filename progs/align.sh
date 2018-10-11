#!/usr/bin/env bash

trap 'printf "error in line %s\n" ${LINENO}; exit;' ERR

declare -r BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"

declare -r outprefix=/tmp/aligntest
declare -r tmpprefix=${outprefix}_tmp

declare -ar batchfiles=(
  '/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/asdt.bed'
  '/net/p01-c2io-nfs/projects/p33/users/franbe/norment_2018/test/qwe.bed'
)

declare -r opt_mini=1
declare -r opt_refallelesfn=""
declare -r opt_samplewhitelist=""

if [ -f "${outprefix}.bed" -a -f "${outprefix}.bim" -a -f "${outprefix}.fam" ] ; then
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


# out: (rs, chr, cm, bp) 
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
get_variant_info_for_dup_chr_cm_bpa() {
  local -r inputfile="$1"
  sort -k 2 ${inputfile} \
    | uniq -u -f 1 \
    | sort -t ' ' -k 1,1 \
    | join -t ' ' -v1 ${inputfile} -
}


# blacklist of variants
declare -r varblacklist=${tmpprefix}.exclude
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
    plink ${plinkflag} ${batchfiles[$i]/%.bed/.bim} --make-just-bim --out ${extmp}
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
  declare plinkoutputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )
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
  plink $flagformat ${plinkinputfn} ${flagextract} --recode transpose --out ${plinkoutputfn}
  plink $flagformat ${plinkinputfn} ${flagextract} --make-bed         --out ${plinkoutputfn}
  # tab-separate all human-readable plink files
  perl -p -i -e 's/[ \t]+/\t/g' ${plinkoutputfn}.bim
  perl -p -i -e 's/[ \t]+/\t/g' ${plinkoutputfn}.fam
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
  declare plinkinputfn=${tmpprefix}_$( get_unique_filename_from_path ${batchfiles[$i]} )
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
  # get unique chr,cm,bpa entries [sort|uniq]
  # sort on rsnumbers corresponding to unique entries
  # get complementary set of rsnumbers (rsnumbers of duplicated vars) [join -v1 | cut]
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
  # get unique chr,cm,bpa entries [sort|uniq]
  # sort on rsnumbers corresponding to unique entries
  # get complementary set (duplicated vars) [join -v1]
  # get the one with highest cm from each [sort -k 3,3 | sort -u]
  # get their rsnumbers [cut | sort -u]
  get_variant_info_for_dup_chr_cm_bpa ${plinkinputfn}.gpz \
    | sort -t ' ' -k 3,3r \
    | sort -t ' ' -u -k 2,2 -k 4,4 \
    | cut -d ' ' -f 1 \
    | sort -u \
    > ${tmpvardups}
  echo -n "$( wc -l ${tmpvardups} | cut -d ' ' -f 1 ) "
  echo "unique coherent duplicate variants retained in chromosome."
  # get unique chr,cm,bpa entries [sort|uniq]
  # sort on rsnumbers corresponding to unique entries
  # get complementary set (duplicated vars) [join -v1]
  # add rsnumbers of non-retained variants to blacklist [join -v1]
  get_variant_info_for_dup_chr_cm_bpa ${plinkinputfn}.gpz \
    | join -t ' ' -v1 - ${tmpvardups} \
    >> ${batchblacklist}
  rm ${tmpvardups}
  if [ ! -z ${opt_refallelesfn} ] ; then
    [ ! -s ${opt_refallelesfn} ] || exit 1
    echo "matching variants to reference.."
    awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${plinkinputfn}.bim \
      | sort -t $'\t' -k 7,7 \
      | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' ${opt_refallelesfn} - \
      | awk -F $'\t' \
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
              print( $1 ) >>batchblacklist
              total_miss++
            }
            else {
              if ( !gmatchx( $2, $3, $5, $6 ) ) {
                print( $5 ) >>batchblacklist
                total_mism++
              }
              else if ( gflip( $2, $3, $5, $6 ) ) {
                print( $5 ) >>batchfliplist
                total_flip++
              }
            }
          } END{
            print( "total missing: ", total_miss )
            print( "total mismatch: ", total_mism )
            print( "total flipped: ", total_flip )
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
  if [ -s "${batchflipfile}" ] ; then
    plinkflags="${plinkflags} --flip ${batchfliplist}"
  fi
  declare plinkoutputfn=${plinkinputfn}_ready_for_mergeattempt
  # NOTE: if plinkflags are empty we could consider "mv $plinkinputfn $plinkoutputfn"
  plink --bfile ${plinkinputfn} ${plinkflags} --make-bed --out ${plinkoutputfn}
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


exit 0


# init global blacklist
# while true
  # for every batch
    # purge batchfile of global blacklist
  # attempt to merge resulting batchfiles
  # stop on success
  # if blacklist is not empty, tell francesco
  # append to global blacklist (derived from errors) 


while true ; do
  tmpbatchfile=$( mktemp .tmpbatch.XXXXXX )
  echo "created temporary batch file ${tmpbatchfile}."
  echo -n > ${tmpbatchfile}
  for batch in $( cat ${mybatchfile} ) ; do
    echo "processing batch ${batch}.."
    tmpfile=$( mktemp .tmpctrl.XXXXXX )
    echo "created temporary variant control file ${tmpfile}."
    batchctrl=${procdir}/$( basename ${batch} ).fatto
    batchvarblacklist=${procdir}/$( basename ${batch} ).exclude
    batchflipfile=${procdir}/$( basename ${batch} ).flip
    if [[ -s "${varblacklist}" ]] ; then
      sort -u ${varblacklist} >> ${batchvarblacklist}
    fi
    j=-1
    tmpbatch=${batch}
    echo -n "current batch is number "
    for (( i = 0; i < ${#mybatchvec[@]}; i++ )); do
      if [[ "${tmpbatch}" == "${mybatchvec[$i]}" ]] ; then
        j=${i}
      fi
    done
    echo $j
    fflag='--bfile'
    if (( j > 0 )) ; then
      case ${myarchvec[$j]} in
        "${bedhex}" )
          fflag='--bfile' ;;
        "${bcfhex}" )
          fflag='--bcf' ;;
        "${vcfhex}" )
          fflag='--vcf' ;;
      esac
      unset myarchvec[$j]
      unset mybatchvec[$j]
    fi
    if [[ ! -f "${batchctrl}" ]] ; then
      if [[ -s "${tmpbatch}.bim" ]] ; then
        for chr in $( cut -f 1 ${tmpbatch}.bim | sort -u ) ; do
          flagextract="--chr ${chr}"
          # use whitelist if existing
          if [ -f "${varwhitelist}" ] ; then
            flagextract="${flagextract} --extract range ${varwhitelist}"
          fi
          tmpvardups=$( mktemp .tmpchr${chr}.XXXXXX )
          chrbatchprefix=${procdir}/$( basename ${tmpbatch} )_chr${chr}
          while [[ ! -f "${chrbatchprefix}.tped" ]] ; do
            plink $fflag ${tmpbatch} ${flagextract} --recode transpose --out ${chrbatchprefix}
          done
          # extract marker information corresponding to unique position-specific genotype series:
          # identical position-specific genotype series will thereby be ignored as harmless here.
          sort -t ' ' -u -k 4 ${chrbatchprefix}.tped | awk '{
            delete catalog
            for ( k = 5; k <= NF; k++ ) {
              if ( $k > 0 ) catalog[$k] = 1
            }
            asorti( catalog )
            printf( "%s %s %s %s", $1, $2, $3, $4 )
            for ( allele in catalog ) printf( "_%s", catalog[allele] )
            print( "" )
          }' | sort -t ' ' -k 4,4 > ${chrbatchprefix}.gp
          echo "extracting marker information for unique position-specific genotype series.."
          cut -d ' ' -f 3,4 ${chrbatchprefix}.gp | sort -t ' ' | uniq -d | sort -t ' ' -k 2,2 \
          | join -t ' ' -1 4 -2 2 ${chrbatchprefix}.gp - | cut -d ' ' -f 3 | sort -u \
          > ${tmpvardups}
          echo -n "$( wc -l ${tmpvardups} | cut -d ' ' -f 1 ) "
          echo "non-coherent duplicate variants marked for deletion in chromosome ${chr}."
          cat ${tmpvardups} >> ${batchvarblacklist}
          sort -t ' ' -k 2,2 ${chrbatchprefix}.tped | awk '{
            delete catalog
            for ( k = 5; k <= NF; k++ ) {
              if ( $k > 0 ) catalog[$k] = 1
            }
            asorti( catalog )
            printf( "%s %s %s", $2, $3, $4 )
            for ( allele in catalog ) printf( "_%s", catalog[allele] )
            print( "" )
          }' | join -t ' ' -v1 - ${tmpvardups} | sort -t ' ' -k 3,3 > ${chrbatchprefix}.gpz
          echo "extracting coherent variants from the set of duplicate variants.."
          cut -d ' ' -f 3 ${chrbatchprefix}.gpz | uniq -d \
          | join -t ' ' -2 3 - ${chrbatchprefix}.gpz | sort -t ' ' -k 3,3r \
          | sort -t ' ' -u -k 1,1 | cut -d ' ' -f 2 | sort -u > ${tmpvardups}
          echo -n "$( wc -l ${tmpvardups} | cut -d ' ' -f 1 ) "
          echo "unique coherent duplicate variants retained in chromosome ${chr}."
          cut -d ' ' -f 3 ${chrbatchprefix}.gpz | uniq -d \
          | join -t ' ' -2 3 - ${chrbatchprefix}.gpz | cut -d ' ' -f 2 | sort \
          | join -t ' ' -v1 - ${tmpvardups} >> ${batchvarblacklist}
          rm ${tmpvardups}
        done
      fi
      touch ${batchflipfile} ${batchvarblacklist}
      if [[ "${refalleles}" != "" && -s "${tmpbatch}.bim" ]] ; then
        echo "matching variants to reference.."
        awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${tmpbatch}.bim | sort -t $'\t' -k 7,7 \
        | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' ${refalleles} - \
        | awk -F $'\t' $AWKLOCINCLUDE \
        -v batchvarblacklist=${batchvarblacklist} \
        -v batchflipfile=${batchflipfile} \
        --source 'BEGIN{
          total_miss = 0
          total_mism = 0
          total_flip = 0
          printf( "" ) >>batchvarblacklist
          printf( "" ) >>batchflipfile
        } {
          if ( $5 == "-" || $6 == "-" ) {
            print( $1 ) >>batchvarblacklist
            total_miss++
          }
          else {
            if ( !gmatchx( $2, $3, $5, $6 ) ) {
              print( $5 ) >>batchvarblacklist
              total_mism++
            }
            else if ( gflip( $2, $3, $5, $6 ) ) {
              print( $5 ) >>batchflipfile
              total_flip++
            }
          }
        } END{
          print( "total missing: ", total_miss )
          print( "total mismatch: ", total_mism )
          print( "total flipped: ", total_flip )
          close( batchvarblacklist )
          close( batchflipfile )
        }'
      fi
      touch ${batchctrl}
    fi
    plinkflags=${plinkflagsdef}
    sort -u ${batchvarblacklist} > $tmpfile
    mv $tmpfile ${batchvarblacklist}
    sort -u ${batchflipfile} > $tmpfile
    mv $tmpfile ${batchflipfile}
    echo "$( wc -l ${batchvarblacklist} ) variants to be excluded."
    echo "$( wc -l ${batchflipfile} ) variants to be flipped."
    if [[ -s "${batchvarblacklist}" ]] ; then
      plinkflags="${plinkflags} --exclude ${batchvarblacklist}"
    fi
    if [[ -s "${batchflipfile}" ]] ; then
      plinkflags="${plinkflags} --flip ${batchflipfile}"
    fi
    if [[ "${plinkflags}" != "" ]] ; then
      mytag=""
      echo "polishing batch ${batch}.."
      while [[ -f "${tmpbatch}.bed" ]] ; do
        tmpbatch=${procdir}/$( basename ${batch} )${mytag}
        mytag="${mytag}Z"
      done
      plink --bfile ${batch} ${plinkflags} --make-bed --out ${tmpbatch}
      if [[ -f "${tmpbatch}.bim" ]] ; then
        perl -p -i -e 's/[ \t]+/\t/g' ${tmpbatch}.bim
      fi
      if [[ -f "${tmpbatch}.fam" ]] ; then
        perl -p -i -e 's/[ \t]+/\t/g' ${tmpbatch}.fam
      fi
      tmpbatchctrl=${procdir}/$( basename ${tmpbatch} ).fatto
      if [[ -f "${batchctrl}" && "${batchctrl}" != "${tmpbatchctrl}" ]] ; then
        cp ${batchctrl} ${tmpbatchctrl}
      fi
    fi
    if [[ -f "${tmpbatch}.bim" ]] ; then
      perl -p -e 's/[ \t]+/\t/g' ${tmpbatch}.bim | awk -F $'\t' '{
        OFS="\t";
        a[1] = $5; a[2] = $6; asort(a);
        $2 = $1":"$4"_"a[1]"_"a[2];
        print;
      }' ${tmpbatch}.bim > ${tmpfile}
      mv ${tmpfile} ${tmpbatch}.bim
    fi
    if [[ -f "${tmpbatch}.bed" ]] ; then
      myarchvec[$j]=${bedhex}
      mybatchvec[$j]=${tmpbatch}
      echo "${tmpbatch}" >> ${tmpbatchfile}
    fi
  done
  if [[ -s "${tmpbatchfile}" ]] ; then
    mv ${tmpbatchfile} ${mybatchfile}
  elif [[ -f "${tmpbatchfile}" ]] ; then
    rm ${tmpbatchfile}
  fi
  misfile=${outprefix}.missnp
  mischrfile=${outprefix}.mischr
  if (( ${#mybatchvec[*]} > 1 )) ; then
    plink --merge-list ${mybatchfile} --out ${outprefix}
    if [[ -f "${outprefix}.bim" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
    fi
    if [[ -f "${outprefix}.fam" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
    fi
    egrep '^Warning: Multiple' ${outprefix}.log | cut -d ' ' -f 7 | tr -d "'." \
    | sort -u >> ${mischrfile}
    if [[ -s "${mischrfile}" ]] ; then
      sort -u ${mischrfile} >> ${misfile}
    fi
  else
    plink --bfile $( cat ${mybatchfile} ) --make-bed --out ${outprefix}
    if [[ -f "${outprefix}.bim" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.bim
    fi
    if [[ -f "${outprefix}.fam" ]] ; then
      perl -p -i -e 's/[ \t]+/\t/g' ${outprefix}.fam
    fi
  fi
  # are we done?
  if [ ! -f "${misfile}" ] ; then
    break
  fi
  # no! prepare to repeat loop
  mv ${misfile} ${varblacklist}
done

