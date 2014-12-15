#!/bin/bash

list() (
  cd "$1"
  for i in *; do
    if [ -f "$i" ]; then
      echo "$PWD/$i"
    elif [ -d "$i" ]; then
      #echo "dir: $PWD/$i"
      list "$i"
    fi
  done
)

getIncludes () (
	echo $(egrep "include" $1 | grep -oE '["<].*.h.*[">]')
)
	 
allFiles=$(echo $(list $1) | tr ' ' '\n')
for x in $allFiles
do 
	echo "[$x]"
	currentIncludes=$(getIncludes $x)
	for y in $currentIncludes
	do 
		#echo "$y"
		noQuotes=$(echo "$y" | sed -e 's/^"//'  -e 's/"$//')
		echo "$noQuotes"
		
	done
done

