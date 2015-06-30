#!/bin/bash
realpathMac() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

if [ $(uname) == "Linux" ]; then
	exit
fi

#elif [ $(uname) == "Darwin" ]; then
#	echo Darwin
#else 
#	echo $(uname)
#fi
#echo $#
if [ $# -ne 2 ]; then
	echo Need two arguments 
	echo "1) directory with binaries to fix"
	echo "2) location of external directory"
	echo "fixDyLinking_mac.sh [BIN_DIR] [EXTERNAL_DIR]"
exit
fi

binDir=$(realpathMac $1)
externalDir=$(realpathMac $2)

FILES=$binDir/*
for f in $FILES
do
echo $f
#cppcms
install_name_tool -change libcppcms.1.dylib $externalDir/local/cppcms/lib/libcppcms.1.dylib $f

install_name_tool -change libbooster.0.dylib $externalDir/local/cppcms/lib/libbooster.0.dylib $f

#armadillo
install_name_tool -change libarmadillo.4.dylib $externalDir/local/armadillo/lib/libarmadillo.4.dylib $f
install_name_tool -change libarmadillo.5.dylib $externalDir/local/armadillo/lib/libarmadillo.5.dylib $f

#bamtools
install_name_tool -change libbamtools.2.3.0.dylib $externalDir/local/bamtools/lib/bamtools/libbamtools.2.3.0.dylib $f
install_name_tool -change libbamtools.2.4.0.dylib $externalDir/local/bamtools/lib/bamtools/libbamtools.2.4.0.dylib $f

#boost libraries
install_name_tool -change libboost_program_options.dylib $externalDir/local/boost/lib/libboost_program_options.dylib $f
install_name_tool -change libboost_thread.dylib $externalDir/local/boost/lib/libboost_thread.dylib $f
install_name_tool -change libboost_filesystem.dylib $externalDir/local/boost/lib/libboost_filesystem.dylib $f
install_name_tool -change libboost_iostreams.dylib $externalDir/local/boost/lib/libboost_iostreams.dylib $f
install_name_tool -change libboost_regex.dylib $externalDir/local/boost/lib/libboost_regex.dylib $f
install_name_tool -change libboost_serialization.dylib $externalDir/local/boost/lib/libboost_serialization.dylib $f
install_name_tool -change libboost_system.dylib $externalDir/local/boost/lib/libboost_system.dylib $f
done