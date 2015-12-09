#!/bin/bash
REFERENCE=~/projects/reference_genomes/human/hg19.fa

module add bwa

for BAM in $(find ~/projects/reference_datasets/human_sequencing/NA12878/ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/supporting/NA12878 -name '*.bam') ; do
	OUT=$(basename $BAM).sr.bam
	if [ -f $OUT -o -f $OUT.unsorted.bam ] ; then
		echo Skipping $OUT
	else
		echo Aligning $OUT
		BWA_PARAM1=""
		if [[ "$BAM" =~ "moleculo" ]] ; then
			BWA_PARAMS="-x intractg"
		fi
		if [[ "$BAM" =~ "pacbio" ]] ; then
			BWA_PARAMS="-x pacbio"
		fi
		samtools bam2fq $BAM | bwa mem $BWA_PARAMS -t $(nproc) $REFERENCE /dev/stdin | samtools view -b - > $OUT.unsorted.bam && \
		novosort -o $OUT $OUT.unsorted.bam
	fi
done