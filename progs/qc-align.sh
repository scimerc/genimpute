#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

declare  cfg_genomebuild; 
         cfg_genomebuild="$( cfgvar_get genomebuild )"; 
readonly cfg_genomebuild

declare  cfg_refprefix
         cfg_refprefix="$( cfgvar_get refprefix )"
readonly cfg_refprefix

if [ ! -z "${cfg_refprefix}" ] ; then
  declare -a batchfiles=( ${opt_inputfiles} "${cfg_refprefix}.all.haplotypes.bcf.gz" )
else
  declare -a batchfiles=( ${opt_inputfiles} )
fi
readonly batchfiles

# default variant reference file
declare  varreffile=${cfg_refprefix}.all.haplotypes.gpa

declare  opt_refcode
         opt_refcode=$( get_unique_filename_from_path "${cfg_refprefix}.all.haplotypes.bcf.gz" )
readonly opt_refcode

#-------------------------------------------------------------------------------

get_plink_varinfo_blacklist() {
# stdin: plink universal variant info file sorted on chromosome and genomic position; stdout: rs
# spare only one variant from sets of coherent colocalized variants, blacklist all else
  awk '{
    blackflag = 0
    vargp = $1"_"$4
    varid = $2
    split( varid, avec, "_" )
    coloc_with_missing = vargp == curvargp && avec[2] == 0 || avec[3] == 0
    if ( NR > 1 && ( varid == curvarid || coloc_with_missing ) ) {
      missing = 0
      for ( k = 5; k <= NF; k++ ) {
        # if the current entry is not valid increment the missing counter
        # else, if the recorded entry is not valid set it to the current one
        # else, if the recorded entry is different from the current one raise
        # the blackflag
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
        print varid
      }
      # if the current version of varid has fewer missing
      if ( missing < curmissing ) {
        # reset the recorded missing counter for varid
        curmissing = missing
        split( varnames, varlist, SUBSEP )
        # add all names varid went by to blacklist
        for ( n in varlist ) print varlist[n]
      } else print varid # add new varid name to blacklist
      # add new name to varid name string
      varnames = varnames SUBSEP varid
    }
    else {
      # else initialize varid arrays
      curmissing = 0
      varnames = varid
      for ( k = 5; k <= NF; k++ ) {
        variome[k] = $k
        if ( $k <= 0 ) curmissing++
      }
    }
    curvargp = vargp
    curvarid = varid
  }' "${*:-/dev/stdin}"
}

#-------------------------------------------------------------------------------

# input: plink bed and eventual tped file sets
# output: purged and [reference] strand-aligned plink file sets

printf "\
  * For every batch:
    * Create a blacklist with non-coherent variants (tped)
    * Append non reference-compatible variants to blacklist
    * Create a fliplist to align strand-flipped variants to reference
    * Purge batch file of blacklist
    * Align fliplist in batch file
" | printlog 0

if [ -z "${cfg_refprefix}" ] ; then
  declare tmpprefix=${opt_outprefix}_tmp
  varreffile=${tmpprefix}.gpa
  for i in ${!batchfiles[@]} ; do
    declare batchcode=$( get_unique_filename_from_path ${batchfiles[$i]} )
    declare b_inprefix=${opt_inprefix}_batch_${batchcode}
    cut -f 2 ${b_inprefix}.bim \
      | awk -F ':' '{
        OFS="\t"
        split( $2, infovec, "_" )
        print( $1, infovec[1], infovec[2], infovec[3] );
      }'
  done \
    | sort -u | sort -k 1,1n -k 4,4n \
    > ${varreffile}
  unset batchcode
  unset b_inprefix
  unset tmpprefix
fi
for i in ${!batchfiles[@]} ; do
  declare batchcode=$( get_unique_filename_from_path ${batchfiles[$i]} )
  declare b_inprefix=${opt_inprefix}_batch_${batchcode}
  declare b_outprefix=${opt_outprefix}_batch_${batchcode}
  declare tmpprefix=${b_outprefix}_tmp
  declare debuglogfn=${tmpprefix}_debug.log
  # check for hash collisions
  if [ -f "${b_outprefix}.bed" -a -f "${b_outprefix}.bim" -a -f "${b_outprefix}.fam" ]; then
    printf "'%s' already exists. skipping alignment step..\n" "${b_outprefix}.bed"
    printf "try increasing 'numchars' in the hash function if you think this should not happen.\n"
    continue
  fi
  if ls ${tmpprefix}* > /dev/null 2>&1; then
    printf "temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
    exit 1
  fi
  printf "Purge and Align batch ${batchfiles[$i]}\n" | printlog 1
  # define input specific plink settings
  declare tmpvarctrl=${tmpprefix}_varctrl
  declare batchallelemap=${tmpprefix}.allelemap
  declare batchblacklist=${tmpprefix}.blacklist
  declare batchfliplist=${tmpprefix}.fliplist
  declare batchidmap=${tmpprefix}.idmap
  sort -k 1,1 -k 4,4 ${b_inprefix}.tped \
    | get_plink_varinfo_blacklist \
    | sort -u > ${batchblacklist}
  {
    printf "$( wc -l ${batchblacklist} | cut -d ' ' -f 1 ) "
    printf "colocalized variants marked for deletion.\n"
  } | printlog 1
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
  # check usability of reference
  printf "matching variants to reference..\n"
  [ -s "${varreffile}" ] || {
    printf "error: file '%s' is unusable.\n" "${varreffile}" >&2;
    exit 1;
  }
  # get chr:bp strings from bim file and join with the corresponding field of varreffile
  awk -F $'\t' '{ OFS="\t"; $7 = $1":"$4; print; }' ${tmpprefix}_nb.bim \
    | sort -t $'\t' -k 7,7 \
    | join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' <( \
      awk '{
        OFS="\t"
        chr = $1
        if ( chr == "MT" ) chr = 26
        if ( chr == "X" || chr == "XY" ) chr = 23
        if ( chr == "Y" ) chr = 24
        print( chr":"$2, $3, $4 )
      }' ${varreffile} \
      | sort -t $'\t' -k 1,1 \
    ) - | awk -F $'\t' \
      -f ${BASEDIR}/lib/awk/nucleocode.awk \
      -f ${BASEDIR}/lib/awk/genotype.awk \
      -f ${BASEDIR}/lib/awk/gflip.awk \
      -f ${BASEDIR}/lib/awk/gmatch.awk \
      -v batchallelemap=${batchallelemap} \
      -v batchblacklist=${batchblacklist} \
      -v batchfliplist=${batchfliplist} \
      -v batchidmap=${batchidmap} \
      --source 'BEGIN{
          OFS="\t"
          total_miss = 0
          total_mism = 0
          total_flip = 0
          total_ambi = 0
          printf( "" ) >batchallelemap
          printf( "" ) >batchblacklist
          printf( "" ) >batchfliplist
          printf( "" ) >batchidmap
        } {
          if ( $5 == "-" || $6 == "-" ) {
            print( $4 ) >batchblacklist
            blackcatalog[$4] = 1
            total_miss++
          }
          else {
            if ( nucleocode($2) == comp_nucleocode($3) ) {
              print( $4 ) >batchblacklist
              blackcatalog[$4] = 1
              total_ambi++
            }
            else {
              if ( !gmatchx( $2, $3, $5, $6 ) ) {
                print( $4 ) >batchblacklist
                blackcatalog[$4] = 1
                total_mism++
              }
              else {
                if ( gflipx( $2, $3, $5, $6 ) ) {
                  print( $4 ) >batchfliplist
                  flipcatalog[$4] = 1
                  total_flip++
                }
                if ( $4 in idmapcatalog ) {
                  print( $4 ) >batchblacklist
                  delete idmapcatalog[$4]
                  blackcatalog[$4] = 1
                }
                else {
                  newid = $1"_"$5"_"$6
                  if ( newid in idcatalog ) {
                    idcatalog[newid]++
                    newid = newid"_"idcatalog[newid]
                  } else idcatalog[newid] = 1
                  idmapcatalog[$4] = newid
                  if ( $2 == 0 || $3 == 0 ) {
                    if ( $2 == 0 ) nref = 3
                    if ( $3 == 0 ) nref = 2
                    aref = $(nref)
                    if ( $4 in flipcatalog )
                      aref = i_to_A( comp_nucleocode(aref) )
                    if ( $5 == aref ) allelemap[$4] = $6
                    if ( $6 == aref ) allelemap[$4] = $5
                  }
                }
              }
            }
          }
        } END{
          print( "total missing:    ", total_miss )
          print( "total mismatch:   ", total_mism )
          print( "total flipped:    ", total_flip )
          print( "total ambiguous:  ", total_ambi )
          for ( varid in idmapcatalog ) {
            split( varid, avec, "_" )
            if ( avec[2] == 0 ) aref = avec[3]
            if ( avec[3] == 0 ) aref = avec[2]
            print( varid, idmapcatalog[varid] ) >batchidmap
            if ( varid in allelemap )
              print( idmapcatalog[varid], 0, aref, allelemap[varid], aref ) >batchallelemap
          }
          close( batchallelemap )
          close( batchblacklist )
          close( batchfliplist )
          close( batchidmap )
        }' \
    | printlog 1
  # update idmap
  sort -t $'\t' -k 1,1 ${batchidmap} > $tmpvarctrl
  mv $tmpvarctrl ${batchidmap}
  # list unique
  sort -t $'\t' -u ${batchfliplist} > $tmpvarctrl
  mv $tmpvarctrl ${batchfliplist}
  printf "$( wc -l ${batchfliplist} ) variants to be flipped.\n" | printlog 1
  # list unique
  sort -t $'\t' -u ${batchblacklist} | join -v1 -t $'\t' - ${batchidmap} > $tmpvarctrl
  mv $tmpvarctrl ${batchblacklist}
  printf "$( wc -l ${batchblacklist} ) variants to be excluded.\n" | printlog 1

  declare plinkflag=""
  if [ -s "${batchblacklist}" ] ; then
    plinkflag="--exclude ${batchblacklist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  ${plinkexec} \
        --bfile ${tmpprefix}_nb ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_nbb \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nbb.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nbb.fam
  declare plinkflag=""
  if [ -s "${batchfliplist}" ] ; then
    plinkflag="--flip ${batchfliplist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  ${plinkexec} \
        --bfile ${tmpprefix}_nbb ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_nbf \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nbf.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_nbf.fam
  ${plinkexec} \
        --bfile ${tmpprefix}_nbf \
        --update-name ${batchidmap} \
        --make-bed \
        --out ${tmpprefix}_un \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  ${plinkexec} \
        --bfile ${tmpprefix}_un \
        --update-alleles ${batchallelemap} \
        --make-bed \
        --out ${tmpprefix}_una \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  declare plinkflag=''
  parcount=$( awk '$1 == 25' ${tmpprefix}_una.bim | wc -l )
  # split X chromosome variants
  if [ $parcount -eq 0 ] ; then
    plinkflag="--split-x ${cfg_genomebuild} no-fail" 
  fi
  ${plinkexec} \
        --bfile ${tmpprefix}_una ${plinkflag} \
        --make-bed \
        --out ${tmpprefix}_out \
        2>&1 >> ${debuglogfn} \
        | tee -a ${debuglogfn}
  unset plinkflag
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.bim
  sed -i -r 's/[ \t]+/\t/g' ${tmpprefix}_out.fam

  rename ${tmpprefix}_out ${b_outprefix} ${tmpprefix}_out.*

  if [ "${batchcode}" == "${opt_refcode}" ] ; then
    cp ${b_outprefix}.bed ${opt_refprefix}.bed
    cp ${b_outprefix}.bim ${opt_refprefix}.bim
    cp ${b_outprefix}.fam ${opt_refprefix}.fam
  fi

  rm ${tmpprefix}*

  unset batchcode
  unset b_inprefix
  unset b_outprefix
  unset debuglogfn
  unset tmpprefix
done

