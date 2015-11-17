#!/bin/bash
#
# runs caller against bams
#
. common.sh
CALLER=lumpy/0.2.11
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [[ -f ${BAM/.bam/.sr.bam} ]] ; then
		BAM=${BAM/.bam/.sr.bam}
	fi
	if samtools view -H $BAM | grep "ID:bwa" >/dev/null ; then
		if samtools view -H $BAM | grep "ID:bwa" | grep -v '\-M' >/dev/null ; then
			echo "Ok to process" >/dev/null
		else 
			# lumpy requires  need bwa mem without -M
			echo "Skipping $BAM as aligned with the bwa -M flag"
			continue
		fi
	else
		echo "Skipping $BAM as not aligned with bwa (or missing bwa @PG header)"
		continue
	fi
	BLACKLIST=""
	if [[ "$CX_BLACKLIST" != "" ]] ; then
		BLACKLIST="-x $CX_BLACKLIST"
	fi
	CX_BAM=$BAM
	CX_CALLER=$CALLER
	cx_save
	XC_OUTPUT=$CX.vcf
	XC_SCRIPT="module add $CALLER; rm -rf $CX; mkdir $CX 2>/dev/null; cd $CX
	samtools view -b -F 1294 $BAM > $CX/discordants.unsorted.bam
	samtools view -h $BAM \
		| ~/src/$CALLER/scripts/extractSplitReads_BwaMem -i stdin \
		| samtools view -Sb - \
		> $CX/splitters.unsorted.bam
	samtools sort $CX/discordants.unsorted.bam $CX/discordants
	samtools sort $CX/splitters.unsorted.bam $CX/splitters
	lumpyexpress $BLACKLIST -B $BAM -S $CX/splitters.bam -D $CX/discordants.bam -o $XC_OUTPUT
	"
	xc_exec
done

