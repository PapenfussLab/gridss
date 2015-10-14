#!/bin/bash
#
# runs caller against bams
#
. common.sh
CALLER=lumpy/0.2.11
for BAM in $DATA_DIR/*.sc.sr.bam ; do
	cx_load $BAM
	if [[ "$CX_ALIGNER" != bwamem ]] ; then
		echo "Skipping $CX as aligned with $CX_ALIGNER"
		continue
	fi
	if [[ "$CX_ALIGNER_FLAGS" == "-M" ]] ; then
		echo "Skipping $CX as bwa aligned with -M"
		continue
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="module add $CALLER; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	# need bwa mem without -M
	samtools view -b -F 1294 $BAM > $CX/discordants.unsorted.bam
	samtools view -h $BAM \
		| ~/src/$CALLER/scripts/extractSplitReads_BwaMem -i stdin \
		| samtools view -Sb - \
		> $CX/splitters.unsorted.bam
	samtools sort $CX/discordants.unsorted.bam $CX/discordants
	samtools sort $CX/splitters.unsorted.bam $CX/splitters
	lumpyexpress -B $BAM -S $CX/splitters.bam -D $CX/discordants.bam -o $XC_OUTPUT
	"
	xc_exec
done

