#!/bin/bash
# parses the argument list according to the options in the variable $opt_list
declare arg_cnt=0
declare arg_list
valid_opt=""
simple_opt=""
composite_opt=""
if [ "${opt_list}" != "" ]; then
	get_opt_check=$(declare -F | grep "get_opt")
	if [ "${get_opt_check}" != "" ]; then
		# echo "parsing arguments [options: ${opt_list[*]}].." > /dev/stderr
		while [ -n "$*" ]; do
			if [ "${composite_opt}" != "" ]; then
#				echo "  reading composite ${composite_opt}.." > /dev/stderr
				argcomp=`echo "$1" | sed "s/^-${composite_opt}//"`
				# if argcomp isn't empty, what's left can be either garbage or the option's value
				if [ "${argcomp}" != "" ]; then
					carg="${argcomp}"
				else
					shift
					carg="$1"
				fi
				buffer=`echo "${carg}" | { grep -e '^=\+$' || true; }`
				if [ "${buffer}" != "" ]; then
					shift
					carg="$1"
				fi
				get_opt ${composite_opt} "${carg}"
				composite_opt=""
				simple_opt=""
				shift
			elif [ "${simple_opt}" != "" ]; then
#				echo "  reading simple ${simple_opt}.." > /dev/stderr
				get_opt ${simple_opt}
				composite_opt=""
				simple_opt=""
				shift
			else
				# echo "parsing '$*'.." > /dev/stderr
				argopt=`echo "$*" | { grep -e "^-" || true; }`
				# echo "commant line string: '${argopt}'.." > /dev/stderr
				if [ "${argopt}" != "" ]; then
					# reset 'valid_opt' value
					# echo "resetting option flag.." > /dev/stderr
					if [ "${valid_opt}" != "" ]; then
						valid_opt=""
					fi
					# found a flag argument (option)
					for copt in ${opt_list[*]}; do
						# echo "matching '${copt}'.." > /dev/stderr
						# remove the '=' sign from the composite flags
						# echo "purging eventual equal signs.." > /dev/stderr
						eopt=`echo ${copt} | sed "s/=//g"`
						# echo "matching '${eopt}'.." > /dev/stderr
						# look for copt at the beginning of the current argument list substring
						fcheck=`echo "$*" | { grep -e "^-${eopt}" || true; }`
						if [ "${fcheck}" != "" ]; then
							valid_opt=${copt}
							# echo "match found." > /dev/stderr
							if [ "${eopt}" == "${copt}" ]; then
								# copt is a simple flag
								# echo "parsing simple option.." > /dev/stderr
								simple_opt=${eopt}
							else
								# copt is a composite flag
								# echo "parsing composite option.." > /dev/stderr
								composite_opt=${eopt}
							fi
						fi
					done
					if [ "${valid_opt}" == "" ]; then
						arg_cnt=$(( arg_cnt + 1 ))
						arg_list[${arg_cnt}]=$1
						composite_opt=""
						simple_opt=""
						shift
					fi
				else
					arg_cnt=$(( arg_cnt + 1 ))
					arg_list[${arg_cnt}]=$1
					composite_opt=""
					simple_opt=""
					shift
				fi
			fi
		done
	else
		echo " opt_parser.sh: module 'get_opt' not defined: can't parse command line options."
		echo "                please define the 'get_opt' function to define parse instructions."
		exit
	fi
else
	echo " opt_parse: no option list found."
	echo "            please provide it by defining the variable opt_list;"
	echo "            simple options are listed as such, composite options with a trailing '=';"
	echo "            all options should be double quoted and stripped of the leading '-' sign."
	echo "            [e.g.: opt_list=(\"h\" \"verbose\" \"a=\" \"booby_trap=\" ....)]"
	exit
fi
echo -n "${arg_list[*]}"
