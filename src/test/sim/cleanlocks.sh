#!/bin/bash
#
# Remove all lock files
#
# this should only be run if a process has been explicitly killed
#
. common.sh

if ! `which qstat >/dev/null 2>&1` ; then
	echo "Please execute from a cluster head node" 1>&2
	exit 1
fi
for LOCK in $DATA_DIR/*.lock ; do
	if [ "$LOCK" == "$DATA_DIR/*.lock" ] ; then
		echo "No locks present"
		exit
	fi
	ID=`basename $LOCK .lock`
	if [ "0${#ID}" -ne "32" ] ; then
		echo "Unexpected ID $ID: aborting to prevent accidental file deletion"
		exit 1
	fi
	LOCKCONTENT=`cat $LOCK`
	# TODO: a better way of determining if the lock is a cluster job or a hostname
	if qstat "$LOCKCONTENT" </dev/null >/dev/null 2>&1 ; then
		echo "$ID still running on cluster. Terminate job before removing lock"
	elif ssh -o StrictHostKeyChecking=no "$LOCKCONTENT" ps -eF 2>/dev/null | grep "$ID" | grep -v grep ; then
		echo "$ID still running on $LOCKCONTENT. Terminate job before removing lock"
	else
		echo "Removing all files for $ID"
		rm -rf $DATA_DIR/$ID
		rm -rf $DATA_DIR/$ID.*
	fi
done
