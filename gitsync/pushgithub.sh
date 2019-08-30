#!/usr/bin/env bash

if [ "$1" != "" ] ; then
  echo "no public repository specified."
  echo "usage: $( basename $0 ) <path_to_repository>"
  echo " e.g.: $( basename $0 ) https://github.com/scimerc/genimpute.git"
fi

PUBREP="$1"
PRJDIR="$( cd $( dirname $( readlink -f "$0" ) )/../ ; pwd )"

git lfs track ${PRJDIR}/lib/data/genetic_map_b37_withX.txt.gz
git add ${PRJDIR}/lib/data/genetic_map_b37_withX.txt.gz
git add ${PRJDIR}/.gitattributes
git commit -m 'lfs'

git remote add origin_gh_norm ${PUBREP}
git push origin_gh_norm

