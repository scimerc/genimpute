#!/usr/bin/env bash

trap 'printf "===> error in %s line %s\n" $(basename $0) ${LINENO}; exit;' ERR

declare -ra batchfiles=( ${opt_inputfiles} )

declare -r cfg_genomebuild="$( cfgvar_get genomebuild )"
declare -r cfg_refallelesfn="$( cfgvar_get refallelesfn )"

#-------------------------------------------------------------------------------

# stdin: plink tped; stdout: rs
# spare only one variant from sets of coherent colocalized variants, blacklist all else
get_tped_blacklist() {
  awk '{
    varid = $1"_"$4
    if ( (varid,5) in variome ) {
      missing = 0
      for ( k = 5; k <= NF; k++ ) {
        if ( variome[varid,k] <= 0 ) variome[varid,k] = $k
        else if ( $k > 0 ) {
          if ( variome[varid,k] != $k ) {
            variome[varid,k] = "_NA_"
            split( varnames[varid], varlist, SUBSEP )
            for ( varname in varlist ) print varname
            print $2
          }
        }
        if ( $k <= 0 ) missing++
      }
      if ( missing < missome[varid] ) {
        missome[varid] = missing
        split( varnames[varid], varlist, SUBSEP )
        for ( varname in varlist ) print varname
      } else print $2
      varnames[varid] = varnames[varid] SUBSEP $2
    } else {
      missome[varid] = 0
      varnames[varid] = $2
      for ( k = 5; k <= NF; k++ ) {
        variome[varid,k] = $k
        if ( $k <= 0 ) missome[varid]++
      }
    }
  }' "${*:-/dev/stdin}"
}

#-------------------------------------------------------------------------------

# for every batch
  # create a batch-chr-specific blacklist with non-coherent variants
  # append non compatible ref variants to blacklist
  # create temp fliplist
  # purge batchfile of blacklist
  # align fliplist in batchfile
# keep blacklists and fliplists just in case


for i in ${!batchfiles[@]} ; do
  declare batchcode=$( get_unique_filename_from_path ${batchfiles[$i]} )
  declare binprefix=${opt_inprefix}_${batchcode}
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
  echo "purging and aligning batch ${batchfiles[$i]}.."
  # define input specific plink settings
  declare tmpvarctrl=${tmpprefix}_varctrl
  declare batchblacklist=${tmpprefix}.blacklist
  declare batchfliplist=${tmpprefix}.fliplist
  get_tped_blacklist ${binprefix}.tped | sort -u > ${batchblacklist} 
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
    awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${binprefix}.bim \
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

  declare plinkflag=""
  if [ -s "${batchblacklist}" ] ; then
    plinkflag="--exclude ${batchblacklist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  plink --bfile ${binprefix} ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_nb \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nb.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nb.fam
  declare plinkflag=""
  if [ -s "${batchfliplist}" ] ; then
    plinkflag="--flip ${batchfliplist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  plink --bfile ${tmpprefix}_nb ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_nbf \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nbf.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nbf.fam
  # rename variants to universal code chr:bp_a1_a2
  make_variant_names_universal_in_bim_file ${tmpprefix}_nbf.bim
  # pre-process sex chromosomes variants
  declare plinkflag=''
  parcount=$( awk '$1 == 25' ${tmpprefix}_nbf.bim | wc -l )
  if [ $parcount -eq 0 ] ; then
    plinkflag="--split-x ${cfg_genomebuild} no-fail" 
  fi
  plink --bfile ${tmpprefix}_nbf ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_out \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
  mv ${tmpprefix}_out.bed ${boutprefix}.bed
  mv ${tmpprefix}_out.bim ${boutprefix}.bim
  mv ${tmpprefix}_out.fam ${boutprefix}.fam
  unset batchcode
  unset binprefix
  unset boutprefix
  unset tmpprefix
  unset debuglogfn
done

rm ${tmpprefix}*

