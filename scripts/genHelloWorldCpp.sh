#!/bin/bash

outFile=""
if [ "$#" -ne 1 ]; then
     outFile="hello.cpp"
else
	outFile=$1
fi

echo "#include <iostream>" >> $outFile
echo "#include <string>" >> $outFile
echo "#include <unistd.h>" >> $outFile
echo "#include <vector>" >> $outFile
echo "#include <stdint.h>" >> $outFile
echo "#include <stdio.h>" >> $outFile
echo "#include <cstddef>" >> $outFile
echo "#include <utility>" >> $outFile
echo "" >> $outFile
echo "int main(int argc, char* argv[])" >> $outFile
echo "{" >> $outFile
echo "	std::cout << \"Hello World!\" << std::endl;" >> $outFile
echo "	return 0;" >> $outFile
echo "}" >> $outFile
