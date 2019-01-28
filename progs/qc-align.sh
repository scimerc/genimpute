#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare  cfg_genomebuild; 
         cfg_genomebuild="$( cfgvar_get genomebuild )"; 
readonly cfg_genomebuild

declare  cfg_refprefix
         cfg_refprefix="$( cfgvar_get refprefix )"
readonly cfg_refprefix

declare -ra batchfiles=( ${opt_inputfiles} "${cfg_refprefix}.all.haplotypes.bcf.gz" )
declare -r  refcode=$( get_unique_filename_from_path "${cfg_refprefix}.all.haplotypes.bcf.gz" )

#-------------------------------------------------------------------------------

# stdin: plink tped sorted on chromosome and genomic position; stdout: rs
# spare only one variant from sets of coherent colocalized variants, blacklist all else
get_tped_blacklist() {
  awk '{
    blackflag = 0
    varid = $1"_"$4
    if ( NR > 1 && varid == curvarid ) {
      missing = 0
      for ( k = 5; k <= NF; k++ ) {
        # if the current entry is not valid increment the missing counter
        # else, if the recorded entry is not valid set it to the current one
        # else, if the recorded entry is different from the current one
        if ( $k <= 0 ) missing++
        else if ( variome[k] <= 0 ) variome[k] = $k
        else if ( variome[k] != $k ) {
          variome[k] = "_NA_"
          blackflag = 1
        }
      }
      if ( blackflag == 1 ) {
        split( varnames, varlist, SUBSEP )
        # add current list of names varid goes by to blacklist
        for ( n in varlist ) print varlist[n]
        # add new name to blacklist
        print $2
      }
      # if the current version of varid has fewer missing
      if ( missing < curmissing ) {
        # reset the recorded missing counter for varid
        curmissing = missing
        split( varnames, varlist, SUBSEP )
        # add all names varid went by to blacklist
        for ( n in varlist ) print varlist[n]
      } else print $2 # add new varid name to blacklist
      # add new name to varid name string
      varnames = varnames SUBSEP $2
    }
    else {
      # else initialize varid arrays
      curmissing = 0
      varnames = $2
      for ( k = 5; k <= NF; k++ ) {
        variome[k] = $k
        if ( $k <= 0 ) curmissing++
      }
    }
    curvarid = varid
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
  declare b_inprefix=${opt_inprefix}_batch_${batchcode}
  declare b_outprefix=${opt_outprefix}_batch_${batchcode}
  declare tmpprefix=${b_outprefix}_tmp
  declare debuglogfn=${tmpprefix}_debug.log
  # check for hash collisions
  if [ -f "${b_outprefix}.bed" -a -f "${b_outprefix}.bim" -a -f "${b_outprefix}.fam" ]; then
    printf "'%s' already exists - skipping recode step..\n" "${b_outprefix}.bed"
    printf "try increasing 'numchars' in the hash function if you think this should not happen.\n"
    continue
  fi
  if ls ${tmpprefix}* > /dev/null 2>&1; then
    printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
    exit 1
  fi
  echo "purging and aligning batch ${batchfiles[$i]}.."
  # define input specific plink settings
  declare tmpvarctrl=${tmpprefix}_varctrl
  declare batchblacklist=${tmpprefix}.blacklist
  declare batchfliplist=${tmpprefix}.fliplist
  sort -k 1,1 -k 4,4 ${b_inprefix}.tped | get_tped_blacklist | sort -u > ${batchblacklist} 
  echo -n "$( wc -l ${batchblacklist} | cut -d ' ' -f 1 ) "
  echo "colocalized variants marked for deletion."
  # use ref alleles if specified or if we have more than one batch
  if [ -z "${cfg_refprefix}" ] ; then
    printf "error: reference is not set.\n" >&2;
    exit 1
  fi
  echo "matching variants to reference.."
  [ -s "${cfg_refprefix}.all.haplotypes.gpa" ] || {
    printf "error: file '%s' is unusable.\n" "${cfg_refprefix}.all.haplotypes.gpa" >&2;
    exit 1;
  }
  # get chr:bp strings from bim file and join with the corresponding field of refprefix
  awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${b_inprefix}.bim \
    | sort -t $'\t' -k 7,7 \
    | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' <( \
      awk '{ OFS="\t"; print( $1":"$2, $3, $4 ); }' ${cfg_refprefix}.all.haplotypes.gpa \
      | sort -t $'\t' -k 1,1 - \
    ) - | awk -F $'\t' \
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
  # list unique
  sort -u ${batchblacklist} > $tmpvarctrl
  mv $tmpvarctrl ${batchblacklist}
  echo "$( wc -l ${batchblacklist} ) variants to be excluded."

  declare plinkflag=""
  if [ -s "${batchblacklist}" ] ; then
    plinkflag="--exclude ${batchblacklist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  ${plinkexec} \
        --bfile ${b_inprefix} ${plinkflag} \
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
  ${plinkexec} \
        --bfile ${tmpprefix}_nb ${plinkflag} \
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
  ${plinkexec} \
        --bfile ${tmpprefix}_nbf ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_out \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam
  mv ${tmpprefix}_out.bed ${b_outprefix}.bed
  mv ${tmpprefix}_out.bim ${b_outprefix}.bim
  mv ${tmpprefix}_out.fam ${b_outprefix}.fam
  if [ "${batchcode}" == "${refcode}" ] ; then
    mv ${b_outprefix}.bed ${opt_refprefix}.bed
    mv ${b_outprefix}.bim ${opt_refprefix}.bim
    mv ${b_outprefix}.fam ${opt_refprefix}.fam
  fi
  rm ${tmpprefix}*
  unset batchcode
  unset b_inprefix
  unset b_outprefix
  unset debuglogfn
done

