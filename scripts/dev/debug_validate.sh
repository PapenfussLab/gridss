#!/bin/bash
for f in $(find . -name '*.sv.bam') ; do
	echo "Seaching for repeated sam records in $f"
	~/dev/gridss/scripts/dev/print_repeats.sh $f
done
echo "Sanity checking evidence"
~/dev/gridss/scripts/gridss.sh \
	-j ~/dev/gridss/target/gridss-*-gridss-jar-with-dependencies.jar \
	-a placeholder.assembly.bam \
	-r masked_ref.fa \
	-w . \
	-o placeholder.gridss.vcf \
	-s preprocess \
	--sanityCheck \
	*.bam