#!/bin/bash
. common.sh

for FILE in $(ls_reference_vcf) ; do
	echo $FILE is a reference VCF
done
for FILE in $(ls_aligned_bam) ; do
	echo $FILE is a sorted bam
done
for FILE in $(ls_aligned_bam "sq") ; do
	echo $FILE is a query sorted bam
done
for FILE in $(ls_aligned_bam "sc" "bwasw") ; do
	echo $FILE is a coordinate sorted bwasw bam
done