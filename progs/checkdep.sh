# check dependencies

{
  locale
  echo
  awk --version 2> /dev/null || { echo 'awk is not installed. aborting..'; exit 1; }
  echo
  join --version 2> /dev/null || { echo 'join is not installed. aborting..'; exit 1; }
  echo
} | printlog 3

