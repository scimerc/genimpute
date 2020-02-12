#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

source "${BASEDIR}/progs/checkdep.sh"

declare -a batchfiles=( ${opt_inputfiles} )
readonly batchfiles

if [ ! -s "${opt_biofile}" ] ; then
  declare N=0 # total number of individuals
  echo -e "> initializing sample biography file..\n"
  # initialize sample biography file
  for i in ${!batchfiles[@]} ; do
    declare batchcode=$( get_unique_filename_from_path "${batchfiles[$i]}" )
    declare b_famfile="${opt_inprefix}_batch_${batchcode}.fam"
    declare Ntmp=$( cat ${b_famfile} | wc -l )
    cut -f 1-4 "${b_famfile}" \
      | awk -F $'\t' \
        -v uid=${cfg_uid} \
        -v batch="${batchcode}" \
        -f "${BASEDIR}/lib/awk/idclean.awk" \
        --source '
          BEGIN{
            OFS="\t";
            print( uid, "BATCH", "FID", "IID", "P1ID", "P2ID" )
          } {
            uid=idclean( $1"_"$2 )
            gsub( "_+$", "", uid )
            gsub( "^_+", "", uid )
            print( uid, batch, $0 )
          }
        '
    N=$(( N + Ntmp ))
    unset Ntmp
    unset b_famfile
    unset batchcode
  done | sort -u -k 1,1 > "${opt_biofile}"
  Nuid=$( cat "${opt_biofile}" | wc -l )
  if [ ${N} -ge ${Nuid} ] ; then
    echo '========> conflicts generated in universal IDs.'
    echo '> please recode IDs so they do not lead to conflicts.'
    exit 1
  fi
fi

