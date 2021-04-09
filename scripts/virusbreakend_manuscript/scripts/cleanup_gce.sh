#!/bin/bash
while true ; do
	echo "$(date) Awake"
	tokill=$( gcloud compute instances list --filter="name~virusbreakend" | grep TERMINATED | cut -f 1 -d ' ' )
	if [[ "$tokill" != "" ]] ; then
		gcloud compute instances delete -q $( gcloud compute instances list --filter="name~virusbreakend" | grep TERMINATED | cut -f 1 -d ' ' )
	else
		echo "Nothing to kill"
	fi
	sleep 180
done