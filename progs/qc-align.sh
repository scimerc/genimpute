#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare  cfg_infosep
         cfg_infosep="$( cfgvar_get infosep )"
readonly cfg_infosep

if [ ! -z "${cfg_refprefix}" ] ; then
  declare -a batchfiles=( ${opt_inputfiles} "${cfg_refprefix}.all.bcf.gz" )
else
  declare -a batchfiles=( ${opt_inputfiles} )
fi
readonly batchfiles

# default variant reference file
declare  varreffile="${cfg_refprefix}.all.gpa"

declare  opt_refcode
         opt_refcode=$( get_unique_filename_from_path "${cfg_refprefix}.all.bcf.gz" )
readonly opt_refcode

#-------------------------------------------------------------------------------

get_plink_varinfo_blacklist() {
# stdin: plink universal variant info file sorted on chromosome and genomic position; stdout: rs
# spare only one variant from sets of coherent colocalized variants, blacklist all else
  awk -v infosep="${cfg_infosep}" '{
    blackflag = 0
    vargp = $1"_"$4
    varid = $2
    split( varid, avec, infosep )
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

echo -e "==== Align ====\n" | printlog 0

printf "\
  * For every batch:
    * Create a blacklist with non-coherent co-localized variants (tped)
    * Append non reference-compatible variants to blacklist
    * Create a fliplist to align strand-flipped variants to reference
    * Purge batch file of blacklist
    * Align fliplist in batch file
\n" | printlog 0

if [ -z "${cfg_refprefix}" ] ; then
  varreffile="${opt_outprefix}.gpa"
  [ -s "${varreffile}" ] || {
    for i in ${!batchfiles[@]} ; do
      declare batchcode=$( get_unique_filename_from_path "${batchfiles[$i]}" )
      declare b_inprefix="${opt_inprefix}_batch_${batchcode}"
      cut -f 2 "${b_inprefix}.bim" | tr "${cfg_infosep}" "\t"
      unset batchcode
      unset b_inprefix
    done \
      | sort -u -k 1,2 > "${opt_outprefix}_tmp.gpa"
    mv "${opt_outprefix}_tmp.gpa" "${varreffile}"
  }
fi
for i in ${!batchfiles[@]} ; do
  declare batchcode=$( get_unique_filename_from_path "${batchfiles[$i]}" )
  declare b_inprefix="${opt_inprefix}_batch_${batchcode}"
  declare b_outprefix="${opt_outprefix}_batch_${batchcode}"
  declare tmpprefix="${b_outprefix}_tmp"
  # check for hash collisions
  if [ -f "${b_outprefix}.bed" -a -f "${b_outprefix}.bim" -a -f "${b_outprefix}.fam" ]; then
    printf "> '%s' already exists. skipping alignment step..\n" "${b_outprefix}.bed"
    printf "> increase 'numchars' in the hash function if you think this shouldn't happen.\n\n"
    continue
  fi
  if ls "${tmpprefix}"* > /dev/null 2>&1; then
    printf "> temporary files '%s*' found. please remove them before re-run.\n" "${tmpprefix}" >&2
    exit 1
  fi
  # define input specific plink settings
  declare tmpvarctrl="${tmpprefix}_varctrl"
  declare batchallelemap="${tmpprefix}.allelemap"
  declare batchblacklist="${tmpprefix}.blacklist"
  declare batchfliplist="${tmpprefix}.fliplist"
  declare batchidmap="${tmpprefix}.idmap"
  sort -k 1,1 -k 4,4 "${b_inprefix}.tped" \
    | get_plink_varinfo_blacklist \
    | sort -u > "${batchblacklist}"
  {
    printf "> $( wc -l "${batchblacklist}" | cut -d ' ' -f 1 ) "
    printf "colocalized variants marked for deletion.\n"
  } | printlog 1
  declare plinkflag=""
  if [ -s "${batchblacklist}" ] ; then
    plinkflag="--exclude ${batchblacklist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  echo -e "  ${plinkexec##*/} --allow-extra-chr
            --bfile ${b_inprefix} ${plinkflag}
            --make-bed
            --out ${tmpprefix}_nb\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
        --bfile "${b_inprefix}" ${plinkflag} \
        --make-bed \
        --out "${tmpprefix}_nb" \
        2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_nb.bim"
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_nb.fam"
  # check usability of reference
  printf "> matching batch '%s' variants to reference..\n" "$( basename "${batchfiles[$i]}" )"
  [ -s "${varreffile}" ] || {
    printf "> error: file '%s' is unusable.\n" "${varreffile}" >&2;
    exit 1;
  }
  # get chr:bp strings from bim file and join with the corresponding field of varreffile
  join -t $'\t' -a2 -2 7 -o '0 2.5 2.6 2.2 1.2 1.3' -e '-' \
    <( \
      awk -F $'\t' -v infosep="${cfg_infosep}" '{
        OFS="\t"
        gsub( "^[Cc]([Hh][Rr]*)*", "", $1 )
        if ( $1 ~ "[Xx][Yy]" ) $1 = 25
        if ( $1 ~ "^[Xx]"    ) $1 = 23
        if ( $1 ~ "^[Yy]"    ) $1 = 24
        if ( $1 ~ "^[Mm]"    ) $1 = 26
        if ( $1 == "25"      ) $1 = 23
        print $1 infosep $2, $3, $4
      }' "${varreffile}" \
      | sort -t $'\t' -k 1,1 \
    ) \
    <( \
      awk -F $'\t' -v infosep="${cfg_infosep}" '{
        OFS="\t"
        chr = $1
        if ( chr == "25" ) chr = 23
        $7 = chr infosep $4
        print
      }' "${tmpprefix}_nb.bim" \
      | sort -t $'\t' -k 7,7 \
    ) \
    | awk -F $'\t' \
      -f "${BASEDIR}/lib/awk/nucleocode.awk" \
      -f "${BASEDIR}/lib/awk/genotype.awk" \
      -f "${BASEDIR}/lib/awk/gflip.awk" \
      -f "${BASEDIR}/lib/awk/gmatch.awk" \
      -v batchallelemap="${batchallelemap}" \
      -v batchblacklist="${batchblacklist}" \
      -v batchfliplist="${batchfliplist}" \
      -v batchidmap="${batchidmap}" \
      -v infosep="${cfg_infosep}" \
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
          # $1 is the join field, i.e. the genomic position "chr*:bp"
          # $2,$3 are the internal alleles
          # $4 is the internal variant name (should be chr:bp:a0:a1)
          # $5,$6 are the reference alleles
          #
          # are any reference alleles missing?
          if ( $5 == "-" || $6 == "-" || $5 <= 0 || $6 <= 0 ) {
            print( $4 ) >batchblacklist
            total_miss++
          }
          else {
            # mark self-complementary (i.e. ambiguous) variants for removal
            if ( nucleocode($2) == comp_nucleocode($3) ) {
              print( $4 ) >batchblacklist
              total_ambi++
            }
            else {
              # mark mismatching colocalized variants (also genuine multi-allelic ones!) for
              # removal. if a matching variant exists in the reference, this will be added to the
              # catalog of valid variants, which takes precedence over the blacklist. note however
              # that this only works if multi-allelic variants are encoded as multiple bi-allelic
              # variants in the reference (which need not be the case for vcf files).(*)
              if ( !gmatchx( $2, $3, $5, $6 ) ) {
                print( $4 ) >batchblacklist
                total_mism++
              }
              else {
                # mark strand-complementary variants for flip
                if ( gflipx( $2, $3, $5, $6 ) ) {
                  print( $4 ) >batchfliplist
                  flipcatalog[$4] = 1
                  total_flip++
                }
                # delete any homonymous variants and mark them for removal
                if ( $4 in idmapcatalog || $4 in idremovecatalog ) {
                  print( $4 ) >batchblacklist
                  idremovecatalog[$4] = 1
                  delete idmapcatalog[$4]
                }
                else {
                  # add variant to catalog
                  # (*) variants in the catalog will be preserved even if later marked for removal
                  newid = $1 infosep $5 infosep $6
                  if ( newid in idcatalog ) {
                    idcatalog[newid]++
                    newid = newid infosep idcatalog[newid]
                  } else idcatalog[newid] = 1
                  idmapcatalog[$4] = newid
                  # if any alleles are missing, use the reference
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
          print( "> total ambiguous:    ", total_ambi )
          print( "> total flipped:      ", total_flip )
          print( "> total mismatch(*):  ", total_mism )
          print( "> total missing:      ", total_miss )
          print( "> ---------------------------------------------------------------------------" )
          print( "> (*) if genuine multi-allelic variants among these are encoded as multiple" )
          print( ">     bi-allelic variants in the reference, they should be preserved." )
          for ( varid in idmapcatalog ) {
            split( varid, avec, infosep )
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
  sort -t $'\t' -k 1,1 "${batchidmap}" > "${tmpvarctrl}"
  mv "${tmpvarctrl}" "${batchidmap}"
  # list unique variants to be flipped
  sort -t $'\t' -u "${batchfliplist}" > "${tmpvarctrl}"
  mv "${tmpvarctrl}" "${batchfliplist}"
  printf "> $( wc -l "${batchfliplist}" ) variants to be flipped.\n" | printlog 1
  # list unique variants to be excluded after cross-cheking against idmap
  sort -t $'\t' -u "${batchblacklist}" | join -v1 -t $'\t' - "${batchidmap}" > "${tmpvarctrl}"
  mv "${tmpvarctrl}" "${batchblacklist}"
  printf "> $( wc -l "${batchblacklist}" ) variants to be excluded.\n" | printlog 1

  declare plinkflag=""
  if [ -s "${batchblacklist}" ] ; then
    plinkflag="--exclude ${batchblacklist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  echo -e "  ${plinkexec##*/} --allow-extra-chr
            --bfile ${tmpprefix}_nb ${plinkflag}
            --make-bed
            --out ${tmpprefix}_nbb\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
        --bfile "${tmpprefix}_nb" ${plinkflag} \
        --make-bed \
        --out "${tmpprefix}_nbb" \
        2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_nbb.bim"
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_nbb.fam"
  declare plinkflag=""
  if [ -s "${batchfliplist}" ] ; then
    plinkflag="--flip ${batchfliplist}"
  fi
  # NOTE: if plinkflags are empty we could consider "mv $opt_batchinpprefix $opt_batchoutprefix"
  echo -e "  ${plinkexec##*/} --allow-extra-chr
            --bfile ${tmpprefix}_nbb ${plinkflag}
            --make-bed
            --out ${tmpprefix}_nbf\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
        --bfile "${tmpprefix}_nbb" ${plinkflag} \
        --make-bed \
        --out "${tmpprefix}_nbf" \
        2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  unset plinkflag
  # tab-separate all human-readable plink files
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_nbf.bim"
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_nbf.fam"
  echo -e "  ${plinkexec##*/} --allow-extra-chr
            --bfile ${tmpprefix}_nbf
            --update-name ${batchidmap}
            --make-bed
            --out ${tmpprefix}_un\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
        --bfile "${tmpprefix}_nbf" \
        --update-name "${batchidmap}" \
        --make-bed \
        --out "${tmpprefix}_un" \
        2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  echo -e "  ${plinkexec##*/} --allow-extra-chr
            --bfile ${tmpprefix}_un
            --update-alleles ${batchallelemap}
            --make-bed
            --out ${tmpprefix}_una\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
        --bfile "${tmpprefix}_un" \
        --update-alleles "${batchallelemap}" \
        --make-bed \
        --out "${tmpprefix}_una" \
        2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  echo -e "  ${plinkexec##*/} --allow-extra-chr
            --bfile ${tmpprefix}_una
            --freq
            --out ${tmpprefix}_una\n" | printlog 2
  ${plinkexec} --allow-extra-chr \
        --bfile "${tmpprefix}_una" \
        --freq \
        --out "${tmpprefix}_una" \
        2> >( tee "${tmpprefix}.err" ) | printlog 3
  if [ $? -ne 0 ] ; then
    cat "${tmpprefix}.err"
  fi
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_una.bim"
  sed -i -r 's/[ \t]+/\t/g' "${tmpprefix}_una.fam"

  rename "${tmpprefix}_una" "${b_outprefix}" "${tmpprefix}_una".*

  if [ "${batchcode}" == "${opt_refcode}" ] ; then
    cp "${b_outprefix}.bed" "${opt_refprefix}.bed"
    cp "${b_outprefix}.bim" "${opt_refprefix}.bim"
    cp "${b_outprefix}.fam" "${opt_refprefix}.fam"
    cp "${b_outprefix}.frq" "${opt_refprefix}.frq"
  fi

  rm "${tmpprefix}"*

  unset batchcode
  unset b_inprefix
  unset b_outprefix
  unset tmpprefix
done

