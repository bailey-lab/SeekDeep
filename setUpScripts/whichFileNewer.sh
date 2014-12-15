#!/bin/bash

fileMod1=$(gstat -c %Y $1)
fileMod2=$(gstat -c %Y $2)

if [ "$fileMod1" -lt "$fileMod2" ]; then
	echo "$2 is newer than $1"
else 
	echo "$1 is newer than $2"
fi