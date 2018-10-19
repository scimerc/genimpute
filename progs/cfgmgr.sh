export _CFGVAR_PREFIX_NAME='_CFGVAR_NAME_'
export _CFGVAR_PREFIX_STATUS='_CFGVAR_STATUS_'
export _CFGVAR_STATUS_EMPTY=0
export _CFGVAR_STATUS_READONLY=1


cfgvar_is_defined() {
  local -r name=$1
  if [ -z "${name}" ] ; then
    printf 'empty variable name\n' >&2
    return 1
  fi
  [ ! -z "${!name+x}" ]
}
export -f cfgvar_is_defined

cfgvar_init() {
  local -r name=$1
  local -r value="$2"
  if cfgvar_is_defined ${_CFGVAR_PREFIX_NAME}$name; then
    printf 'variable %s is already defined.\n' $name >&2
    return 1
  fi
  export ${_CFGVAR_PREFIX_NAME}${name}="${value}"
  export ${_CFGVAR_PREFIX_STATUS}${name}=${_CFGVAR_STATUS_EMPTY}
}
export -f cfgvar_init

_cfgvar_set_status() {
  local -r name=$1
  local -r status=$2
  local -r varname=${_CFGVAR_PREFIX_NAME}${name}
  if ! cfgvar_is_defined $varname; then
    printf 'variable %s is not defined.\n' $name >&2
    return 1
  fi
  local -r statusname=${_CFGVAR_PREFIX_STATUS}${name}
  if [ -z "${statusname}" ] ; then
    printf 'undefined variable status.\n' >&2
    return 1
  fi
  export ${_CFGVAR_PREFIX_STATUS}${name}=${status}
}
export -f _cfgvar_set_status


cfgvar_set_readonly() {
  _cfgvar_set_status "$1" $_CFGVAR_STATUS_READONLY
}
export -f cfgvar_set_readonly

_cfgvar_read_file() {
  local -r filename=$1
  local -r func=$2
  while read line; do
    #TODO: include tabs in the following
    # remove leading spaces
    line="$( echo "${line}" | sed "s/^\ *//" )"
    # remove trailing spaces
    line="$( echo "${line}" | sed "s/\ *$//" )"
    # skip if line is empty
    if [ -z "${line}" ] ; then
      continue
    fi
    # skip if line starts with
    if [[ "${line}" =~ ^# ]] ; then
      continue
    fi
    echo "$line"
    # split line in name and value
    name=${line%%=*}
    value=${line#*=}
    # do stuff with name and value
    ${func} "${name}" "${value}"
  done < <( cat ${filename} )
}
export -f _cfgvar_read_file

cfgvar_init_from_file() {
  _cfgvar_read_file $1 cfgvar_init
}
export -f cfgvar_init_from_file

cfgvar_update_from_file() {
  _cfgvar_read_file $1 cfgvar_set
}
export -f cfgvar_update_from_file

_cfgvar_list_names() {
  env | grep "^${_CFGVAR_PREFIX_NAME}" | sed "s/^${_CFGVAR_PREFIX_NAME}//" | sed 's/=.*$//'
}
export -f _cfgvar_list_names

cfgvar_setall_readonly() {
  for name in $( _cfgvar_list_names ) ; do
    cfgvar_set_readonly ${name}
  done
}
export -f cfgvar_setall_readonly

cfgvar_get() {
  local -r usrname=$1
  local -r varname=${_CFGVAR_PREFIX_NAME}$1
  if ! cfgvar_is_defined $varname; then
    printf 'variable %s is not defined.\n' $usrname >&2
    return 1
  fi
  echo "${!varname}"
}
export -f cfgvar_get

cfgvar_set() {
  local -r name=$1
  local -r value="$2"
  local -r varname=${_CFGVAR_PREFIX_NAME}$name
  if ! cfgvar_is_defined $varname; then
    printf 'variable %s is not defined.\n' $name >&2
    return 1
  fi
  local -r statusname=${_CFGVAR_PREFIX_STATUS}${name}
  if [ "${!statusname}" -eq $_CFGVAR_STATUS_READONLY ] ; then
    printf 'variable %s is read-only.\n' $name >&2
    return 1
  fi
  export ${varname}="${value}"
}
export -f cfgvar_set

cfgvar_show_config() {
  todo
}
export -f cfgvar_show_config

