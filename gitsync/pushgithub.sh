#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

if [ "$*" == "" ] ; then
  echo "no public repository specified."
  echo "usage: $( basename $0 ) <path_to_repository>"
  echo " e.g.: $( basename $0 ) https://github.com/scimerc/genimpute.git"
  exit 0
fi

PUBREP="$1"
PRJDIR="$( cd $( dirname $( readlink -f "$0" ) )/../ ; pwd )"

cd ${PRJDIR}

echo 'committing changes..'

git add .
git commit -m 'repository export'

echo 'adding large files support..'

git lfs track lib/data/genetic_map_b37_withX.txt.gz
git add lib/data/genetic_map_b37_withX.txt.gz
git add .gitattributes
git commit -m 'lfs' || true

echo 'pushing changes to remote..'

git remote add origin_gh_norm ${PUBREP}
git pull origin_gh_norm master
git commit -m 'merge with remote'
git push origin_gh_norm

