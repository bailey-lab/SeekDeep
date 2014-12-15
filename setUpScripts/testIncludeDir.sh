#!/bin/bash
outFile=$2
echo "digraph G  { " > $outFile
echo "bgcolor =\"#000000\" ">> $outFile
echo "labelloc=\"t\" " >> $outFile
echo "fontcolor = \"#ffffff\"" >> $outFile
echo "fontsize = 20 " >> $outFile
echo "label = \"top\" ">> $outFile
echo "fixedsize = true; " >> $outFile
for x in $(gfind $1 -iregex '.*\(\.c\|\.cpp\|\.h\|\.hpp\)'); do
	if [[ $(basename $x) == *.h*  ]] 
	then
		 echo $(echo $(basename $x) | gsed 's/\./_/g') [shape=circle,style=filled,fixesize =true, color = \"#000000\", fillcolor =\"#0571b0\", width = 1]>>$outFile ;
	else 
		echo $(echo $(basename $x) | gsed 's/\./_/g') [shape=circle,style=filled,fixesize =true, color = \"#000000\", fillcolor =\"#ca0020\", width = 1]>>$outFile ;
	fi
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


