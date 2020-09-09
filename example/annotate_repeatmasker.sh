#!/bin/bash

# gridss48f_4b is a good example of where the detailed RepeatMasker alignment
# provides information not present in the summary.
# For most purposes, the .out file is sufficient and the detailed alignments
# are not necessary.
../scripts/gridss_annotate_vcf_repeatmasker.sh \
	-j ../target/gridss-2.10.0-gridss-jar-with-dependencies.jar \
	-w out \
	-o out/annotate_repeatmasker_example.annotated.vcf \
	../src/test/resources/repeatmasker/ebv/merged.bwa.bam.gridss.vcf \
	--rmargs "-species human -a"