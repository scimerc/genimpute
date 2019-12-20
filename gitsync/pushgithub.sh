#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

declare -r MAXFILESIZE=50 #MB
declare -r REMOTEREPOS='origin_gh_norm'

if [ "$*" == "" ] ; then
  echo "no public repository specified."
  echo "usage: $( basename $0 ) <path_to_repository>"
  echo " e.g.: $( basename $0 ) https://github.com/scimerc/genimpute.git"
  exit 0
fi

PUBREP="$1"
PRJDIR="$( cd $( dirname $( readlink -f "$0" ) )/../ ; pwd )"

cd ${PRJDIR}

lflist=$( find . -type f -size +${MAXFILESIZE}M -not -path '*.git*'  )

if [ "${lflist}" != "" ] ; then

  echo 'adding large files support..'

  git lfs track ${lflist}
  git add .gitattributes ${lflist}
  git commit -m 'lfs' || true

fi

echo 'committing and pushing changes to remote..'

git add .
git remote -v | grep -qw "${REMOTEREPOS}" || git remote add ${REMOTEREPOS} ${PUBREP}
git pull ${REMOTEREPOS} master
git commit -m 'merge with remote'
git push ${REMOTEREPOS}

