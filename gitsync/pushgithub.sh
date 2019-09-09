#!/usr/bin/env bash

if [ "$1" == "" ] ; then
  echo "no public repository specified."
  echo "usage: $( basename $0 ) <path_to_repository>"
  echo " e.g.: $( basename $0 ) https://github.com/scimerc/genimpute.git"
  exit 0
fi

PUBREP="$1"
PRJDIR="$( cd $( dirname $( readlink -f "$0" ) )/../ ; pwd )"

cd ${PRJDIR}

git lfs track lib/data/genetic_map_b37_withX.txt.gz
git add lib/data/genetic_map_b37_withX.txt.gz
git add .gitattributes
git commit -m 'lfs'

git remote add origin_gh_norm ${PUBREP}
git push origin_gh_norm

