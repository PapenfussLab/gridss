#!/bin/bash
getopt --test
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
	echo 'WARNING: "getopt --test"` failed in this environment.' 1>&2
	echo "WARNING: The version of getopt(1) installed on this system might not be compatible with the GRIDSS driver script." 1>&2
fi
set -o errexit -o pipefail -o noclobber -o nounset
last_command=""
current_command=""
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "\"${last_command}\" command completed with exit code $?."' EXIT
#253 forcing C locale for everything
export LC_ALL=C

EX_USAGE=64
EX_NOINPUT=66
EX_CANTCREAT=73
EX_CONFIG=78

dbname=virusbreakenddb
USAGE_MESSAGE="
VIRUSBreakend database build

Downloads and builds the kraken2 database used by VIRUSBreakend

This program requires kraken2 and samtools to be on PATH

Usage: virusbreakend-build.sh [--db $virusbreakenddb]"
OPTIONS=h
LONGOPTS=help,db:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit $EX_USAGE
fi
eval set -- "$PARSED"
while true; do
	case "$1" in
		-h|--help)
			echo "$USAGE_MESSAGE" 1>&2
			exit 0
			;;
		--db)
			dbname=virusbreakenddb
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 1
			;;
	esac
done
for tool in kraken2-build samtools gunzip wget awk tar dustmasker rsync ; do
	if ! which $tool >/dev/null; then
		echo "Error: unable to find $tool on \$PATH" 1>&2
		exit $EX_CONFIG
	fi
	write_status "Found $(which $tool)"
done

kraken2-build --download-taxonomy --db $dbname
kraken2-build --download-library human --db $dbname
kraken2-build --download-library viral --db $dbname
kraken2-build --download-library UniVec_Core --db $dbname

cd $dbname
wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.genomic.fna.gz
# download virushostdb files
# convert virushostdb.genome.fna to kraken2 notation
gunzip -c virushostdb.genomic.fna.gz | awk -f <(awk '
BEGIN { FS = "\t"; print "BEGIN {\n" }
{ split($4, contig, ", "); for (i in contig) { print "\tlookup[\"[" contig[i] "]\"] = \">kraken:taxid|" $1 "|" "\" ;" } }
END { print "}\n{ if (substr($1, 1, 1) == \">\" && lookup[$2] != \"\") { $1=lookup[$2]substr($1, 2) } ; print }" } ' < virushostdb.tsv) > virusbreakend.virushostdb.genomic.fna
cd -
kraken2-build --add-to-library $dbname/virusbreakend.virushostdb.genomic.fna --db $dbname
# TODO why does masking result in empty .mask files?
kraken2-build --build --db $dbname
for f in $(find $dbname/ -name '*.fna') ; do samtools faidx $f; done

tar -czvf virusbreakend.db.$dbname.tar.gz \
	$dbname/*.k2d \
	$dbname/taxonomy/nodes.dmp \
	$dbname/virushostdb.tsv \
	$dbname/library/viral/*.fna* \
	$dbname/library/added/*.fna* \

