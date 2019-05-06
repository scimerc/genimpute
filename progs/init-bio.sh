#!/usr/bin/env bash

# exit on error
set -Eeou pipefail

source "${BASEDIR}/progs/checkdep.sh"

if [ ! -s "${opt_biofile}" ] ; then
  echo "> initializing sample biography file.."
  # initialize sample biography file
  cut -f 1-4 "${opt_outprefix}.fam" \
    | awk -F $'\t' -v uid=${cfg_uid} -f "${BASEDIR}/lib/awk/idclean.awk" \
      --source '
        BEGIN{
          OFS="\t";
          print( uid, "FID", "IID", "P1ID", "P2ID" )
        } { print( idclean( $1"_"$2 ), $0 ) } \
      ' \
    | sort -u -k 1,1 > "${opt_biofile}"
  Norg=$( cat "${opt_outprefix}.fam" | wc -l )
  Nuid=$( cat "${opt_biofile}" | wc -l )
  if [ ${Norg} -ge ${Nuid} ] ; then
    echo '========> conflicts generated in universal IDs.'
    echo '> please recode IDs so they do not lead to conflicts.'
    exit 1
  fi
fi

