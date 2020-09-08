#!/bin/bash
GRIDSS_JAR=$PWD/../target/gridss-2.10.0-gridss-jar-with-dependencies.jar
input_vcf=$PWD/../src/test/resources/repeatmasker/ebv/merged.bwa.bam.gridss.vcf
mkdir -p out && \
java -Xmx64m -cp $GRIDSS_JAR \
	gridss.InsertedSequencesToFasta \
	I= $input_vcf \
	O= out/annotate_repeatmasker_example.fa && \
RepeatMasker -a -pa $(nproc) -species human -dir out out/annotate_repeatmasker_example.fa && \
# We can annotate with full alignment information by using the .align file
# GRIDSS can parse both the .fa.out and the .fa.align files
java -Xmx1g -cp $GRIDSS_JAR \
	gridss.repeatmasker.AnnotateVariantsRepeatMasker \
	I= $input_vcf \
	O= out/annotate_repeatmasker_example.annotated.vcf \
	RM= out/annotate_repeatmasker_example.fa.align
