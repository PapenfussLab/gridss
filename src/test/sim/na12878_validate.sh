#!/bin/bash
SUPPORTING=~/projects/reference_datasets/human_sequencing/NA12878/ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/supporting/NA12878/

SUPPORTING_SR=~/na12878/longread

for F in $(ls $SUPPORTING_SR/*.sr.bam) ; do
	java -Xmx16g \
		-cp ~/bin/gridss-0.8.1-SNAPSHOT-jar-with-dependencies.jar \
		au.edu.wehi.validation.SplitReadToBed \
		I=$F \
		O=/dev/stdout \
	| bedtools sort > $F.bed
done


#java -Xmx16g \
	-cp ~/bin/gridss-0.8.1-SNAPSHOT-jar-with-dependencies.jar \
	au.edu.wehi.validation.ValidateDeletions \
	SCM=216 \
	SWS=216 \
	I=data.na12878/tovalidate.bedpe \
	O=data.na12878/tovalidate.annotated.bedpe \
	LR=$SUPPORTING/moleculo/alignment/NA12878.moleculo.bwa-mem.20140110.bam \
	LR=$SUPPORTING/pacbio/alignment/chemistry1.sorted.bam \
	LR=$SUPPORTING/pacbio/alignment/chemistry_2_picard.bam \
	LR=$SUPPORTING/pacbio/alignment/chemistry_3_picard.bam \
	LR=$SUPPORTING/pacbio/alignment/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam
