#!/bin/bash
jvm_args=" \
		-Dpicard.useLegacyParser=false \
		-Dsamjdk.use_async_io_read_samtools=true \
		-Dsamjdk.use_async_io_write_samtools=true \
		-Dsamjdk.use_async_io_write_tribble=true \
		-Dsamjdk.buffer_size=4194304 \
		-Dsamjdk.async_io_read_threads=$(nproc)"
input_bam=~/dev/colo829/COLO829R_dedup.realigned.bam
for n in 10000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000; do
	if [[ ! -f n_${n}_.insert_size_metrics ]] ; then
		echo Process \"$n\"
		java -Xmx2g $jvm_args \
			-cp ~/dev/gridss/target/gridss-2.10.0-gridss-jar-with-dependencies.jar \
			gridss.analysis.CollectGridssMetrics \
			--INPUT $input_bam \
			--OUTPUT "n_${n}_" \
			--THRESHOLD_COVERAGE 50000 \
			--TMP_DIR . \
			--FILE_EXTENSION null \
			--STOP_AFTER $n
	fi
done

grep -A 100000 "## HISTOGRAM" *.insert_size_metrics | tr '_-' '\t' | cut -d '	' -f 2,6,7 | grep -E "	[0-9]+	" > insertsize.tsv
grep "S$" *.cigar_metrics | tr '_:\-' '\t' | cut -d '	' -f 2,5,6 > softclip.tsv
grep -A 1 MAX_READ_LENGTH *.idsv_metrics | grep idsv_metrics | tr ':-' '\t' | sed 's/^n_//'| sed 's/_.idsv_metrics//' | sed 's/.*MAX_READ_LENGTH/n\tMAX_READ_LENGTH/' | sort -r | uniq > idsv.tsv
