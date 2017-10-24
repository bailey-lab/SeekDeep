#!/usr/bin/env bash

if [ $# -ne 1 ]; then
	echo "Need to supply number of CPUs to use, e.g. ./install.sh 7"
	exit
fi

#make sure on the master branch and then pull for any new commits
git checkout master
git pull
#re-run install
./configure.py
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk --overWrite
make -j $1
