#!/bin/bash
checkPackage(){
    retval=""
    output=$(dpkg -l | grep -E '^ii' | cut -d' ' -f3 | grep '^'"$1"'$')
    if [ ! -z "$output" -a "$output" != " " ]; then
	retval="true"
	#echo "$(tput setaf 2)$(tput bold)Package $1 is installed$(tput sgr0)"
    else 
	retval="false"
	#echo "$(tput setaf 1)$(tput bold)Package $1 is not installed$(tput sgr0)"
    fi
    echo "$retval"
}

pkgInstalled=$(checkPackage $1)

if [ "$pkgInstalled" == "true" ]
then
     echo "$(tput setaf 2)$(tput bold)Package $1 is installed$(tput sgr0)"
else
     echo "$(tput setaf 1)$(tput bold)Package $1 is not installed$(tput sgr0)"
fi 

