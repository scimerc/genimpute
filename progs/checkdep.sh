# check dependencies

{
  locale
  echo
  awk  --version 2> /dev/null || { echo 'awk is not installed. aborting..' ; exit 1; }
  echo
  gzip --version 2> /dev/null || { echo 'gzip is not installed. aborting..'; exit 1; }
  echo
  sed  --version 2> /dev/null || { echo 'sed is not installed. aborting..' ; exit 1; }
  echo
} | printlog 1

