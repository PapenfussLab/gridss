#!/bin/bash
dir=publicdata/sim
for caller in gridss manta svaba novobreak ; do
	outfile=$dir/merged_${caller}_sim.vcf
	grep -E "^##" $dir/gen/run1_/run1_.$caller.vcf > $outfile
	echo '##INFO=<ID=RUNID,Number=1,Type=String,Description="simulation run">' >> $outfile
	grep "#CHROM" $dir/gen/run1_/run1_.$caller.vcf >> $outfile
	grep -vE "^#" $(find $dir -name "*.$caller.vcf.full.vcf") | awk -v OFS="\t" '{ split($1, s, ":"); $1=s[2]; $8 = "RUNID=" s[1] ";" $8 ; print }' >> $outfile
done
# for gridss:
#grep -vE "^#" $(find $dir -name "*.$caller.vcf.full.vcf") | awk -v OFS="\t" '{ split($1, s, ":"); $1=s[2]; $8 = "RUNID=" s[1] ";" $8 ; print }' >> $outfile