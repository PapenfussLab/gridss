#!/bin/bash
# 
# Moves all results for the given caller to an archive subdirectory
#

if [[ "$1" == "" ]] ; then
	echo "Caller arg missing"
	exit
fi

echo "Archiving $1"

for DATA in data.* ; do
#for DATA in data.model ; do
	echo "echo Processing $DATA"
	cd $DATA
	mkdir -p archive
	for F in $(grep "$1" *.metadata) ; do
		if [[ "$2" != "go" ]] ; then
			echo "$F"
			echo "mv ${F:0:32}* archive"
		else
			mv ${F:0:32}* archive
		fi
	done
	cd ..
done

