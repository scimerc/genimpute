#!/usr/bin/env bash
# exports repository to a gzip-compressed tarball

# exit on error
set -ETeuo pipefail

[ -d ".git" ] || { printf "error at %d.\n" $LINENO >&2; exit 1; }

[ "$*" != "" ] || { printf "usage: %s <dir>\n" ${BASH_SOURCE[0]}; exit 0; }

EXPDIR="$1"

LOCTMPDIR=/tmp

PRJNAME='genimpute'
PRJDIR="$( cd $( dirname $( readlink -e "${BASH_SOURCE[0]}" ) )/../ ; pwd )"
PRJDIR=$( readlink -e $PRJDIR )

GITROOT="$(git rev-parse --show-toplevel)"
# calculate md5, but only use first 6 chars
MD5=$(find ${GITROOT}/.git -type f | xargs cat | md5sum)
MD5=${MD5:0:7}

echo "Cloning source repository '${PRJDIR}'..."
cd ${LOCTMPDIR} && rm -rf ${PRJNAME}
git clone ${PRJDIR} ${PRJNAME}

echo "Compressing into '${PRJNAME}.git.tar.gz'..."
tar czf ${PRJNAME}.git.tar.gz ${PRJNAME}

declare -r OUTFILENAME=${PRJNAME}_"$(date +"%y%m%d-%H")"_${MD5}.git.tar.gz
mv ${PRJNAME}.git.tar.gz "${EXPDIR}/${OUTFILENAME}"

echo "Repository exported to '${EXPDIR}/${OUTFILENAME}'."

