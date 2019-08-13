#!/usr/bin/env bash

if [ $# -ne 1 ]; then
	echo "Need to supply number of CPUs to use, e.g. ./upgradeDevelop.sh 7"
	exit
fi

#make sure on the master branch and then pull for any new commits
git checkout develop
git pull
#re-run install
./configure.py
rm -fr external/build/bamtools/ external/build/jsoncpp/ external/build/restbed/ 
rm -fr external/local/bamtools/ external/local/jsoncpp/ external/local/restbed/ 
rm -fr external/local/njhseq/ external/local/seqServer/ external/local/TwoBit/ external/local/njhcpp/

./setup.py --compfile compfile.mk --outMakefile makefile-common.mk --overWrite

make -j $1
