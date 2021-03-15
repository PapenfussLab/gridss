#!/bin/bash
cd ../protecteddata/hmfv4/
grep "^#" $(grep -l BEALN *breakend.vcf | head -1) | grep -v "##contig" | grep -v "#CHROM" > merged.vcf
echo "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample identifier\">" >> merged.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample" >> merged.vcf
grep -vE "^(#|(adjusted))" *.virusbreakend.vcf | sed 's/.virusbreakend.vcf:/	/' | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\tSAMPLE=" $1 ";" $9 "\t" $10 "\t" $11 }' >> merged.vcf
echo -n "sample	" > merged.tsv
cat *.virusbreakend.vcf.summary.tsv | grep taxid_genus | head -1 >> merged.tsv
grep -v taxid_genus *.virusbreakend.vcf.summary.tsv | sed -r 's/^([^.]+).virusbreakend.vcf.summary.tsv:/\1\t/' >> merged.tsv
ls -1 *.vcf | cut -f 1 -d . > samples.tsv
