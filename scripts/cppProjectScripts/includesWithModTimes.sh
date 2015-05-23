#!/bin/bash
outFile=$2
for x in $(gfind $1 -iregex '.*\(\.c\|\.cpp\|\.h\|\.hpp\)'); do
	echo $(echo $(basename $x) | gsed 's/\./_/g') gstat -c Y% $x >> $outFile;
done

for x in $(gfind $1 -iregex '.*\(\.c\|\.cpp\|\.h\|\.hpp\)'); do for y in $(echo $(egrep "#include" $x | grep -oE '["].*\.h.*["]') | sed '/^$/d'| gsed -e 's/\b"//g'  -e 's/"\b//g'); do 
	if [[ $(basename $x) == *.h*  ]] 
	then 
		echo $(echo $(basename $y) \-\> $(basename $x) | gsed 's/\./_/g') [penwidth=5, color=\"#92c5de\"] >>$outFile;
	else
		echo $(echo $(basename $y) \-\> $(basename $x) | gsed 's/\./_/g') [penwidth=5, color=\"#f4a582\"] >>$outFile; 
	fi
	done ;
done

echo "}" >> $outFile


