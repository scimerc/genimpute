#!/usr/bin/env bash

trap 'exit' ERR

[ -d ".git" ] || { printf "error at %d\n" $LINENO >&2; exit 1; }

PRJNAME='tsdtools'
PRJDIR="$(pwd)"

GITROOT="$(git rev-parse --show-toplevel)"
MD5=$(find ${GITROOT}/.git -type f | xargs cat | md5sum)
MD5=${MD5:0:7}

cd /tmp && rm -rf ${PRJNAME}.git
git clone ${PRJDIR} ${PRJNAME}.git

# calculate md5, but only use first 6 chars


printf "compressing...\n"
tar cfz ${PRJNAME}.git.tar.gz ${PRJNAME}.git
printf "done\n"

declare -r OUTFILENAME=${PRJNAME}_"$(date +"%y%m%d-%H")"_${MD5}.git.tar.gz
mv ${PRJNAME}.git.tar.gz ${OUTFILENAME}

printf "\n%s\n" "${OUTFILENAME}"

cp $OUTFILENAME /tsd/p33/fx/export/$(echo $USER | sed "s/^[^-]*-//")/
