#!/bin/bash
#
# Remove all files not associated with a metadata file
#
. common.sh

LIST=""
for FILE in $DATA_DIR/* ; do
	CX=$(try_get_cx_prefix $FILE)
	if [ "$CX" == "" -o ! -f "$CX.metadata" ] ; then
		if [[ -d $FILE ]] ; then
			if [[ $FILE =~ [0-9a-f]{32} ]] ; then
				echo "Removing $FILE" 1>&2
				LIST="$LIST $FILE"
			else
				echo "Skipping directory $FILE"
			fi
		else
			echo "Removing $FILE" 1>&2
			LIST="$LIST $FILE"
		fi
	fi
done
if [ "$LIST" == "" ] ; then
	echo "Nothing to clean" 1>&2
else 
	echo "Sleeping for 5s: last chance to kill this script before they're all gone!" 1>&2
	sleep 5
	rm -r $LIST
fi
