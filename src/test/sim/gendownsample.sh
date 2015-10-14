#!/bin/bash
#
# Generates perfect read data sets
#
. common.sh

# Generate simulated reads for each reference VCF
for BAM in $DATA_DIR/*.sc.bam ; do
	cx_load $BAM
	if [ -z "$CX_DOWNSAMPLE_FROM" ] ; then
		if [[ "$CX_READ_LENGTH" -ne 100 && $FULL_MATRIX -ne 1 ]] ; then
			echo "Skipping read length $CX_READ_LENGTH ($BAM)"
			continue
		fi
		if [[ "$CX_READ_FRAGMENT_LENGTH" -ne 300 && $FULL_MATRIX -ne 1 ]] ; then
			echo "Skipping fragment size $CX_READ_FRAGMENT_LENGTH ($BAM)"
			continue
		fi
		if [[ "$CX_INDEL_REALIGNMENT" != "" ]] ; then
			echo "Skipping indel realigned ($BAM)"
			continue
		fi
		SRC_DEPTH=$CX_READ_DEPTH
		for DEPTH in $READ_DEPTHS ; do
			if [ $DEPTH -lt $SRC_DEPTH ] ; then
				CX_READ_DEPTH=$DEPTH
				CX_DOWNSAMPLE_FROM=$BAM
				cx_save
				PR=$(echo "scale=4; $CX_READ_DEPTH / $SRC_DEPTH" | bc) 
				XC_OUTPUT=$CX.sc.bam
				XC_MEMORY=3072
				XC_SCRIPT="
				DownsampleSam I=$BAM O=$CX.sc.downsample.bam PROBABILITY=$PR &&\
				SortSam I=$CX.sc.downsample.bam O=$CX.sq.bam SO=queryname &&\
				mv $CX.sc.downsample.bam $XC_OUTPUT &&\
				samtools index $XC_OUTPUT
				"
				xc_exec
			fi
		done
	fi
done
