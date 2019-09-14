#!/usr/bin/env bash

# exit on error
set -ETeuo pipefail

declare BASEDIR="$( cd "$( dirname $0 )" && cd .. && pwd )"
declare IDIR=$HOME

printf "enter install path [default='%s']: " "${IDIR}"
read -r

[ "${REPLY}" == "" ] || IDIR="${REPLY}"

[ -d "${IDIR}" ] || {
  printf "no such directory: '%s'.\n" "${IDIR}"
  exit 0
}

printf "copying files to '%s'..\n" "${IDIR}"

mkdir -p \
  ${IDIR}/lib/config \
  ${IDIR}/lib/data \
  ${IDIR}/lib/awk \
  ${IDIR}/progs

[ -e "${IDIR}/lib/3rd" ] || { cd "${IDIR}/lib"; ln -s ../progs 3rd; }

cp -ruv ${BASEDIR}/lib/config/* ${IDIR}/lib/config/
cp -ruv ${BASEDIR}/lib/data/* ${IDIR}/lib/data/
cp -ruv ${BASEDIR}/lib/awk/* ${IDIR}/lib/awk/
cp -ruv ${BASEDIR}/lib/3rd/* ${IDIR}/lib/3rd/
cp -ruv ${BASEDIR}/progs/* ${IDIR}/progs/

