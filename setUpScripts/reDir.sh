#!/bin/bash

#for x in `gfind . -type d`; do echo -e "${x}" ; done

recreateDir()(
	for x in `find $1 -type d`; do echo -e $(echo "${x}" | sed "s/${1//\//\/}/${2//\//\/}/g" ) ; done	
)

recreateDir $1 $2

for x in `gfind $1 -iregex '.*\(\.h\|\.hpp\)'`; do echo -e $(echo "${x}" | sed "s/${1//\//\/}/${2//\//\/}/g" ) ; done 

