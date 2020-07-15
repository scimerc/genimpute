#!/usr/bin/env bash
# imports repository files from a tarball

# exit on error
set -ETeuo pipefail

[ -d ".git" ] || { printf "error at %d.\n" $LINENO >&2; exit 1; }

[ "$*" != "" ] || { printf "usage: %s <tarfile>\n" ${BASH_SOURCE[0]}; exit 0; }

TARFILE="$1"

LOCTMPDIR=/tmp

PRJNAME='genimpute'
PRJDIR="$( cd $( dirname $( readlink -e "${BASH_SOURCE[0]}" ) )/../ ; pwd )"
PRJDIR=$( readlink -e $PRJDIR )

SELFABS=$( readlink -e "${BASH_SOURCE[0]}" )
SELFREL="${SELFABS#${PRJDIR}/}"

echo "Extracting source archive '${TARFILE}'..."
cd ${LOCTMPDIR} && rm -rf ${PRJNAME}
tar xvf "${TARFILE}"
echo "done."

[ -d ${PRJNAME} ] || { printf "project '%s' not found.\n" ${PRJNAME}; exit 1; }

echo "Syncing into '${PRJDIR}'..."
rsync -av --delete --exclude={".git","${SELFREL}"} ${PRJNAME}/ ${PRJDIR}/
echo "done."

[ -f ${PRJNAME}/${SELFREL} ] || { echo "${SELFREL}: Not daring to self-destroy."; exit 0; }
cmp $SELFABS ${PRJNAME}/${SELFREL} 2>/dev/null || { echo "${SELFREL}: Not daring to self-sync."; }

