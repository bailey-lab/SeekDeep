#!/usr/bin/env bash

if [ $# -ne 1 ]; then
	echo "Need to supply number of CPUs to use, e.g. ./install.sh 7"
	exit
fi

./configure.py
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk 
make -j 7
