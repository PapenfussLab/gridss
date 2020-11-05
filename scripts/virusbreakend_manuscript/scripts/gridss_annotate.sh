#!/bin/bash
export PATH=~/dev/gridss/scripts:$PATH
export GRIDSS_JAR=~/dev/gridss/target/gridss-2.10.1-gridss-jar-with-dependencies.jar

hcc_dir=../publicdata/hcc/
mkdir -p $hcc_dir/gridss/
for f in $hcc_dir/gridss*.vcf ; do
	gridss_annotate_vcf_kraken2.sh \
		-o $hcc_dir/gridss/$(basename $f).k2.vcf \
		-t 12 \
		--kraken2db ~/dev/virusbreakend/virusbreakenddb \
		$f
	grep INSTAXID $hcc_dir/gridss/$(basename $f).k2.vcf | grep -v "INSTAXID=9606;" > $hcc_dir/gridss/$(basename $f).viral.vcf
done


