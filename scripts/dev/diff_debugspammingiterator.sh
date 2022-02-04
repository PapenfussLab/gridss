#!/bin/bash

for f1 in DebugSpammingIterator$1.Worker*.log ; do
	f2=${f1/tor$1/tor$2}
	# let truncated files still be happy
	lines1=$(cat $f1 | wc -l)
	lines2=$(cat $f2 | wc -l)
	lines=$lines2
	if [[ $lines1 -lt $lines2 ]] ; then
		lines=$lines1
	fi
	for filter in \
			"SAMEvidenceSource.iterator().rawiterator" \
			"SAMEvidenceSource.shouldFilterPreTransform" \
			"SAMEvidenceSource.transform" \
			"SAMEvidenceSource.iterator(QueryInterval[])" \
			"SAMEvidenceSource.shouldFilter(SAMRecord)" \
			"SAMEvidenceSource.shouldFilter(DirectedEvidence)" \
			"asEvidence" \
			"SAMEvidenceSource.iterator(QueryInterval[]).filter" \
			"mergedIterator" \
			"PositionalAssembler.evidenceIt" \
			"PositionalAssembler.SupportNodeIterator" \
			"PositionalAssembler.AggregateNodeIterator" \
			"PositionalAssembler.PathNodeIterator" \
			"AssemblyEvidenceSource.assembler" \
			; do
		(diff --speed-large-files -q \
			<(head -$lines $f1 | grep -F "$filter") \
			<(head -$lines $f2 | grep -F "$filter") >/dev/null \
			&& echo "---------------IDENTICAL	$filter	$f1") || \
			   echo "----------NOT--IDENTICAL	$filter	$f1"
		diff --speed-large-files\
			<(head -$lines $f1 | grep -F "$filter") \
			<(head -$lines $f2 | grep -F "$filter") | head -2
	done
done


