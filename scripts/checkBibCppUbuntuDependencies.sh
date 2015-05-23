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
declare -a arr=("git" "wget" "subversion" "cmake" "g++-4.8" "gcc-4.8" "libbz2-dev" "python-dev" "libopenmpi-dev" "r-base-dev" "libcurl4-openssl-dev")
notInstalled=()
echo "Need the following"
echo ${arr[*]}
COUNTER=0
for i in "${arr[@]}"
do
    #echo "$COUNTER"
    pkgInstalled=$(checkPackage $i)
    if [ "$pkgInstalled" == "true" ]
    then
      	echo "$(tput setaf 2)$(tput bold)Package $i is installed$(tput sgr0)"
    else
     	echo "$(tput setaf 1)$(tput bold)Package $i is not installed$(tput sgr0)"
	notInstalled+=($i)
    fi 
    COUNTER=$[$COUNTER +1]
done
tLen=${#notInstalled[@]}
if [ $(echo $tLen) -ne "0" ]
then
    echo "Need these packages"
    echo ${notInstalled[*]}
    echo "do"
    echo "sudo apt-get install ${notInstalled[*]}"
else
    echo "All set for dependencies now run setUp.py while in bib-cpp directory" 
fi


