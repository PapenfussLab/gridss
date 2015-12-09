#!/bin/bash

for F in chemistry1.sorted.bam.sr.bam chemistry_2_picard.bam.sr.bam chemistry_3_picard.bam.sr.bam NA12878.moleculo.bwa-mem.20140110.bam.sr.bam NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sr.bam ; do
	if [[ ! -f $F.sr.bed ]] ; then
		java -ea -Xmx16g -cp ~/bin/gridss-*-jar-with-dependencies.jar au.edu.wehi.validation.BamToBed \
			I=$F \
			SR=$F.sr.bed \
			SP=null \
			OD=true \
			MIN_FINAL_INDEL_SIZE=25 &
	fi
done

for F in chemistry1.sorted.bam chemistry_2_picard.bam chemistry_3_picard.bam NA12878.moleculo.bwa-mem.20140110.bam  NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam ; do
	if [[ ! -f $F.sp.bed ]] ; then
		java -ea -Xmx16g -cp ~/bin/gridss-*-jar-with-dependencies.jar au.edu.wehi.validation.BamToBed \
			I=$F \
			SR=null \
			SP=$F.sp.bed \
			OD=true \
			MIN_FINAL_INDEL_SIZE=25 &
	fi
done
