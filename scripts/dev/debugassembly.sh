#!/bin/bash
# Creates a minimal GRIDSS working set for debugging a managable subset of the data
#
# $1 source working directory
# $2 directory to create
# $3 file contains reads to extract. One per line
function usage {
	echo "debugchunk.sh <sourceDir>  <destDir> <chr:start-end>"
	exit 1
}
if [[ ! -d $1 ]] ; then
	usage
fi
if [[ "$2" == "" ]] ; then
	usage
fi
if [[ "$3" == "" ]] ; then
	usage
fi
SRC=$1
DEST=$2
FILTER=$3
cd $SRC
for W in *.gridss.working ; do
	cd $SRC/$W
	mkdir -p $DEST/$W
	echo Processing $W
	BAM=$(basename $W .gridss.working)
	cp $SRC/$W/$BAM.alignment_summary_metrics $DEST/$W
	cp $SRC/$W/$BAM.cigar_metrics $DEST/$W
	cp $SRC/$W/$BAM.coverage.blacklist.bed $DEST/$W
	cp $SRC/$W/$BAM.idsv_metrics $DEST/$W
	cp $SRC/$W/$BAM.insert_size_metrics $DEST/$W
	cp $SRC/$W/$BAM.mapq_metrics $DEST/$W
	cp $SRC/$W/$BAM.quality_distribution_metrics $DEST/$W
	cp $SRC/$W/$BAM.sv_metrics $DEST/$W
	cp $SRC/$W/$BAM.tag_metrics $DEST/$W
	samtools view -b $SRC/$W/$BAM.sv.bam $FILTER > $DEST/$W/$BAM.sv.bam
	samtools index $DEST/$W/$BAM.sv.bam
	ln -s $DEST/$W/$BAM.sv.bam $DEST/$BAM
done
