export _CFGVAR_PREFIX_NAME_='_CFGVAR_NAME_'
export _CFGVAR_PREFIX_STATUS_='_CFGVAR_STATUS_' # stores read-only flag
export _CFGVAR_PREFIX_RCOUNT_='_CFGVAR_RCOUNT_' # read counter
export _CFGVAR_STATUS_EMPTY_=0
export _CFGVAR_STATUS_READONLY_=1

#-------------------------------------------------------------------------------

# get variable value
cfgvar_get() {
  local -r name=$1
  _cfgvar_get ${name}
  if [ $? -ne 0 ] ; then
    return 1
  fi
  # increase read counter
  local -r rcvar=${_CFGVAR_PREFIX_RCOUNT_}${name}
  export $rcvar=$((rcvar + 1))
}
export -f cfgvar_get

#-------------------------------------------------------------------------------

# initialize variable
cfgvar_init() {
  local -r name=$1
  local -r value="$2"
  if cfgvar_is_defined ${_CFGVAR_PREFIX_NAME_}$name; then
    printf 'variable %s is already defined.\n' $name >&2
    return 1
  fi
  export ${_CFGVAR_PREFIX_NAME_}${name}="${value}"
  export ${_CFGVAR_PREFIX_STATUS_}${name}=${_CFGVAR_STATUS_EMPTY_}
  export ${_CFGVAR_PREFIX_RCOUNT_}${name}=0
}
export -f cfgvar_init

#-------------------------------------------------------------------------------

# initialize variable values from configuration file
cfgvar_init_from_file() {
  _cfgvar_read_file $1 cfgvar_init
}
export -f cfgvar_init_from_file

#-------------------------------------------------------------------------------

# check if variable is defined
cfgvar_is_defined() {
  local -r name=$1
  if [ -z "${name}" ] ; then
    printf 'empty variable name\n' >&2
    return 1
  fi
  [ ! -z "${!name+x}" ]
}
export -f cfgvar_is_defined

#-------------------------------------------------------------------------------

# set variable value
cfgvar_set() {
  local -r name=$1
  local -r value="$2"
  local -r varname=${_CFGVAR_PREFIX_NAME_}$name
  if ! cfgvar_is_defined $varname; then
    printf 'variable %s is not defined.\n' $name >&2
    return 1
  fi
  local -r statusname=${_CFGVAR_PREFIX_STATUS_}${name}
  if [ "${!statusname}" -eq $_CFGVAR_STATUS_READONLY_ ] ; then
    printf 'variable %s is read-only.\n' $name >&2
    return 1
  fi
  export ${varname}="${value}"
}
export -f cfgvar_set

#-------------------------------------------------------------------------------

# set variable to read-only
cfgvar_set_readonly() {
  _cfgvar_set_status "$1" $_CFGVAR_STATUS_READONLY_
}
export -f cfgvar_set_readonly

#-------------------------------------------------------------------------------

# set all variables to read-only
cfgvar_setall_readonly() {
  for name in $( _cfgvar_list_names ) ; do
    cfgvar_set_readonly ${name}
    if [ $? -ne 0 ] ; then
      return 1
    fi
  done
}
export -f cfgvar_setall_readonly

#-------------------------------------------------------------------------------

# print all variable and their values
cfgvar_show_config() {
  for name in $( _cfgvar_list_names ); do
    # use internal get() so we do not increaase read counter
    printf "%s = %s\n" ${name} "$(_cfgvar_get ${name})"
  done
}
export -f cfgvar_show_config

#-------------------------------------------------------------------------------

# read variable values from configuration file
cfgvar_update_from_file() {
  _cfgvar_read_file $1 cfgvar_set
}
export -f cfgvar_update_from_file

#-------------------------------------------------------------------------------

# undefine all variables and functions
cfgvar_cleanup() {
  # check read counters XXX: disabled
  #for name in $( _cfgvar_list_names ) ; do
  #  rcvar=${_CFGVAR_PREFIX_RCOUNT_}${name}
  #  if [ ${!rcvar} -eq 0 ] ; then
  #    printf "warning: unused variable '%s'\n" ${name} >&2
  #  fi
  #done
  # cleanup vars
  for name in $(_cfgvar_list_varnames_prefix "_CFGVAR_") ; do
    unset ${name}
  done
  # get all function names with _cfgvar or cfgvar prefix
  local -r funcs="$(typeset -F | cut -d " " -f 3 | grep "^cfgvar_\|^_cfgvar_")"
  for name in $funcs; do
    unset ${name}
  done
}
export -f cfgvar_cleanup

#-------------------------------------------------------------------------------

_cfgvar_get() {
  local -r usrname=$1
  local -r varname=${_CFGVAR_PREFIX_NAME_}${usrname}
  if ! cfgvar_is_defined $varname; then
    printf 'variable %s is not defined.\n' $usrname >&2
    return 1
  fi
  printf "%s" "${!varname}"
}
export -f _cfgvar_get

#-------------------------------------------------------------------------------

_cfgvar_read_file() {
  local -r filename="$1"
  local -r function=$2
  if [ ! -f "$filename" ]; then
    printf "error: file '%s' not found\n" ${filename} >&2
    return 1
  fi
  local linecounter=0
  while read line; do
    linecounter=$((linecounter + 1))
    #TODO: include tabs in the following
    # remove leading spaces
    line="$( echo "${line}" | sed 's/^\ *//g;' )"
    # remove trailing spaces
    line="$( echo "${line}" | sed 's/\ *$//g;' )"
    # skip if line is empty
    if [ -z "${line}" ] ; then
      continue
    fi
    # skip if line starts with
    if [[ "${line}" =~ ^# ]] ; then
      continue
    fi
    # TODO if line does not contain a '=' return an error
    # remove spaces at end of name and at beginning of value
    line="$( printf '%s' "${line}" | sed 's/\ *=\ */=/g;' )"
    # split line in name and value
    name=${line%%=*}
    # TODO if name is empty return an error
    value=${line#*=}
    # do stuff with name and value
    value="$( printf '%s' "${value}" | sed 's/^["'\'']//g;''s/["'\'']$//g;' )"
    ${function} "${name}" "${value}"
    if [ $? -ne 0 ] ; then
      return 1
    fi
  done < <( cat ${filename} )
}
export -f _cfgvar_read_file

#-------------------------------------------------------------------------------

_cfgvar_set_status() {
  local -r name=$1
  local -r status=$2
  local -r varname=${_CFGVAR_PREFIX_NAME_}${name}
  if ! cfgvar_is_defined $varname; then
    printf 'variable %s is not defined.\n' $name >&2
    return 1
  fi
  local -r statusname=${_CFGVAR_PREFIX_STATUS_}${name}
  if [ -z "${statusname}" ] ; then
    printf 'undefined variable status.\n' >&2
    return 1
  fi
  export ${_CFGVAR_PREFIX_STATUS_}${name}=${status}
}
export -f _cfgvar_set_status

#-------------------------------------------------------------------------------

_cfgvar_list_varnames_prefix() {
  local -r pfx="$1"
  env | grep "^${pfx}" | sed 's/=.*$//'
}
export -f _cfgvar_list_varnames_prefix

#-------------------------------------------------------------------------------

_cfgvar_list_names() {
  _cfgvar_list_varnames_prefix ${_CFGVAR_PREFIX_NAME_} \
    | sed "s/^${_CFGVAR_PREFIX_NAME_}//"
}
export -f _cfgvar_list_names

