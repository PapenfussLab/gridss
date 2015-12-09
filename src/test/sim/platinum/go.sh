#!/bin/bash
OUTPUT=platinum-na12878.vcf
REFERENCE=~/projects/reference_genomes/human/hg19.fa

#CleanSam I=~/projects/reference_datasets/human_sequencing/NA12878/ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12878_S1.bam O=ERA172924_NA12878_S1.bam
#CleanSam I=~/projects/reference_datasets/human_sequencing/NA12878/ftp.sra.ebi.ac.uk/vol1/ERA207/ERA207860/bam/NA12878_S1.bam VALIDATION_STRINGENCY=LENIENT O=/dev/stdout | samtools view -h - | samtools view -b - > tmp.ERA207860_NA12878_S1.bam
#ReorderSam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT I=tmp.ERA207860_NA12878_S1.bam O=ERA207860_NA12878_S1.bam S=true R=$REFERENCE &
#ReorderSam VALIDATION_STRINGENCY=LENIENT I=tmp.ERA207860_NA12878_S1.bam O=/dev/stdout s=true R=~/projects/reference_genomes/human/hg19.fa | samtools view -h - | samtools view -b - | novosort -f -i -o ERA207860_NA12878_S1.sorted.bam -
#ReorderSam VALIDATION_STRINGENCY=LENIENT I=tmp.ERA172924_NA12878_S1.bam O=/dev/stdout s=true R=~/projects/reference_genomes/human/hg19.fa | samtools view -h - | samtools view -b - | novosort -f -i -o ERA172924_NA12878_S1.sorted.bam -

gridss() {
	if [[ -f $OUTPUT ]] ; then
		echo "Already processed"
		exit
	fi
	rm -f realign.sh
	java \
		-ea \
		-Xmx64g \
		-Dgridss.visualisation.savetimeouts=true \
		-Dgridss.visualisation.trackAssemblyProgress=true \
		-jar ~/bin/gridss-*-jar-with-dependencies.jar \
		INPUT=ERA172924_NA12878_S1.sorted.bam \
		TMP_DIR=. \
		WORKING_DIR=. \
		OUTPUT=$OUTPUT \
		REFERENCE=$REFERENCE \
		SCRIPT=realign.sh \
		VERBOSITY=DEBUG \
		MIN_SCORE=100 \
		READ_PAIR_CONCORDANT_PERCENT=0.99 \
		2>&1 | tee -a idsv.std.$$.log
		# mate pair library ERA207860 has Nextera adapters not trimmed
		# INPUT=ERA207860_NA12878_S1.sorted.bam \
		
	if [[ -f realign.sh ]] ; then
		source realign.sh
		gridss
	fi
}
gridss

if [ ! -f rp.bam ] ; then
	echo "Generating debugging read pair bam"
	samtools merge -r rp.bam $(find . | grep .rp.bam) && samtools index rp.bam &
fi
if [ ! -f sc.bam ] ; then
	echo "Generating debugging soft clip bam"
	samtools merge -r sc.bam $(find . | grep .sc.bam) && samtools index sc.bam &
fi
if [ ! -f scr.bam ] ; then
	echo "Generating debugging soft clip realignment bam"
	#samtools merge -r scr.bam $(find . -name '*.realign.bam' | grep .scr.bam) && samtools index scr.bam &
fi
wait
