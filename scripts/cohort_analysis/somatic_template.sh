#!/bin/bash

#
# Example gridss pipeline for somatic tumour/normal analysis
#
BLACKLIST=%(blacklist_file)s
REFERENCE=%(reference_genome)s
OUTPUT=%(patient_id)s.sv.vcf
ASSEMBLY=${OUTPUT/.sv.vcf/.gridss.assembly.bam}
GRIDSS_JAR=%(gridss_jarfile)s

if [[ ! -f $GRIDSS_JAR ]] ; then
	echo "Missing $GRIDSS_JAR. Update the GRIDSS_JAR variable in the shell script to your location"
	exit 1
fi

./gridss.sh \
	--reference $REFERENCE \
	--output $OUTPUT \
	--assembly $ASSEMBLY \
	--jar $GRIDSS_JAR \
	--workingdir . \
	--blacklist $BLACKLIST \
	2>&1 | tee -a gridss.$HOSTNAME.$$.log
