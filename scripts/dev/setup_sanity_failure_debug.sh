#!/bin/bash
# $1: working directory
if [[ ! -d "$1" ]] ; then
	echo "specify working directory"
	exit 1
fi
dir=sanity_failure_debug/$(basename $1)
rm -r $dir 2>/dev/null
mkdir -p $dir
cp $1/*_metrics $dir
cp $1/*.bed $dir
head -1 reads_failing_sanity_check.txt > reads_failing_sanity_check.txt.tmp.txt
picard FilterSamReads I=$(ls -1 $1/*.sv.bam) O=$dir/$(basename $(ls -1 $1/*.sv.bam)) FILTER=includeReadList READ_LIST_FILE=reads_failing_sanity_check.txt
