#!/bin/bash
REFERENCE=~/projects/reference_genomes/human/hg19.fa
SOURCE=~/projects/reference_datasets/human_sequencing/NA12878/ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/
F=ERR194147
if [[ ! -f $F.sorted.bam ]] ; then
	bwa mem -t 20 $REFERENCE $SOURCE/${F}_1.fastq.gz $SOURCE/${F}_2.fastq.gz | samtools view -b - > $F.unsorted.bam && \
	novosort --md -o $F.sorted.bam $F.unsorted.bam 
fi
if [[ ! -f $F.sc.bam ]] ; then
	AddOrReplaceReadGroups RGPL=illumina RGPU=ERR194147 RGLB=ERR194147 RGCN=ERR194147 RGSM=ERR193147 I=ERR194147.sorted.bam O=ERR194147.sc.bam
fi
# bowtie2 version
if [[ ! -f $F.bt2.sorted.bam ]] ; then
	bowtie2 -X 1000 -p $(nproc) --mm --very-sensitive-local -x $REFERENCE -1 $SOURCE/${F}_1.fastq.gz -2 $SOURCE/${F}_2.fastq.gz | samtools view -b - > $F.bt2.unsorted.bam && \
	novosort --md -o $F.bt2.sorted.bam $F.bt2.unsorted.bam 
fi
if [[ ! -f $F.bt2.sc.bam ]] ; then
	AddOrReplaceReadGroups RGPL=illumina RGPU=ERR193147 RGLB=ERR193147 RGCN=ERR193147 RGSM=ERR193147 I=$F.bt2.sorted.bam O=$F.bt2.sc.bam
fi

#samtools view -h ERA172924_NA12878_S1.sorted.bam | AddOrReplaceReadGroups RGPL=illumina RGPU=platinum RGLB=platinum RGCN=platinum RGSM=platinum I=/dev/stdin O=ERA172924.sc.bam
