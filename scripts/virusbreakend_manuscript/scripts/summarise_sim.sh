#!/bin/bash
cat $(find gen/batvi -name 'final_hits.txt') | grep -v LIB > gen/all_final_hits.txt
grep MSA gen/all_final_hits.txt > ../publicdata/sim/batvi_all_final_hits_MSA.txt
grep -v MSA gen/all_final_hits.txt > ../publicdata/sim/batvi_all_final_hits_noMSA.txt
grep "^#" $(find gen/virusbreakend -name 'virusbreakend_*_*x.vcf' | head -1) > ../publicdata/sim/virusbreakend.vcf
grep -v "##" $(find gen/virusbreakend -name 'virusbreakend_*_*x.vcf' ) | grep -v "#CHROM" >> ../publicdata/sim/virusbreakend.vcf
grep -v Confidence $(find gen/verse -name integration-sites.txt) > ../publicdata/sim/verse.tsv
grep chr $(find gen/vifi -name 'output.clusters.txt.range') | tr ':' ',' > ../publicdata/sim/vifi.tsv
grep "^#" gen/gridss/1/60x/gridss_1_60x.vcf > ../publicdata/sim/gridss2.vcf
grep -v "##" $(find gen/gridss -name 'gridss_*x.vcf') | grep -v "#CHROM" >> ../publicdata/sim/gridss2.vcf
grep chr $(find gen/batvi -name predictions.opt.subopt.txt) | tr ':' '\t' > ../publicdata/sim/batvi_lc.tsv
grep chr $(find gen/batvi/ -name '*.predictionsx.msa.txt') | tr ':' '\t' > ../publicdata/sim/batvi_msa.tsv
grep "^#" $(find gen/gsv -name 'gsv_*_*x.vcf' | head -1) > ../publicdata/sim/gsv.vcf
grep -v "##" $(find gen/gsv -name 'gsv_*x.vcf') | grep -v "#CHROM" | grep -E '\[|\]' >> ../publicdata/sim/gsv.vcf


grep -E "(Start)|(End)" $(find gen -name 'log*out') > ../publicdata/sim/timing.txt
sed 's/^.*\///' ../publicdata/sim/timing.txt | sed 's/.out:/\t/' | tr '._ ' '\t\t\t'| sed 's/\t\t/\t/' | cut -f 2,5,6,10,11,14  > ../publicdata/sim/timing.tsv
