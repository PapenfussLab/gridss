#!/bin/bash
cd ../protecteddata/hmf/
grep "[0-9]" $(find . -name '*viral.txt') | tr ':' '\t' > viral.out
grep "[0-9]" *extracted.txt  | tr ':' '\t' > extracted.out
grep -A 1 GENOME_TERRITORY *wgs_metrics.txt | grep -v GENOME_TERRITORY | grep -v ""[-][-]"" | sed 's/.txt-/	/' > wgs_metrics.out
grep "^#" $(ls -1 *.virusbreakend.vcf | head -1) | grep -v "##contig" | grep -v "#CHROM" > merged.vcf
echo "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample identifier\">" >> merged.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample" >> merged.vcf
grep -vE "^(#|(adjusted))" *.virusbreakend.vcf | sed 's/.virusbreakend.vcf:/	/' | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\tSAMPLE=" $1 ";" $9 "\t" $10 "\t" $11 }' >> merged.vcf