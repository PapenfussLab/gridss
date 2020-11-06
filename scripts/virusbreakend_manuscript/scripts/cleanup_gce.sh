#!/bin/bash
while true ; do
	echo "$(date) Awake"
	gcloud compute instances delete -q $( gcloud compute instances list --filter="name~virusbreakend" | grep TERMINATED | cut -f 1 -d ' ' )
	sleep 180
done
