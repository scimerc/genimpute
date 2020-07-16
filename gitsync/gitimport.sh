#!/usr/bin/env bash
# imports repository files from a tarball

# exit on error
set -ETeuo pipefail

[ -d ".git" ] || { printf "Error at %d.\n" $LINENO >&2; exit 1; }

[ "$*" != "" ] || { printf "Usage: %s <tarfile>\n" ${BASH_SOURCE[0]}; exit 0; }

TARFILE="$1"

LOCTMPDIR=/tmp

PRJNAME='genimpute'
PRJDIR="$( cd $( dirname $( readlink -e "${BASH_SOURCE[0]}" ) )/../ ; pwd )"
PRJDIR=$( readlink -e $PRJDIR )

SELFABS=$( readlink -e "${BASH_SOURCE[0]}" )
SELFREL="${SELFABS#${PRJDIR}/}"

# get list of large files (not to be synced)
LFSLIST=( $( git lfs ls-files | cut -d ' ' -f 3 ) )

echo "Extracting source archive '${TARFILE}'..."
cd ${LOCTMPDIR} && rm -rf ${PRJNAME}
tar xvf "${TARFILE}"

[ -d ${PRJNAME} ] || { printf "Project '%s' not found.\n" ${PRJNAME}; exit 1; }

echo "Syncing into '${PRJDIR}'..."
rsync -av --delete --exclude={".git","${LFSLIST[@]}","${SELFREL}"} ${PRJNAME}/ ${PRJDIR}/

[ -f ${PRJNAME}/${SELFREL} ] || { echo "${SELFREL}: Not daring to self-destroy."; exit 0; }
cmp $SELFABS ${PRJNAME}/${SELFREL} 2>/dev/null || { echo "${SELFREL}: Not daring to self-sync."; }

echo "Repository synced."

