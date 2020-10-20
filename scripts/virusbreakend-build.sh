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

Downloads and builds the databases used by VIRUSBreakend

An the 

This program requires kraken2 and samtools to be on PATH

Usage: virusbreakend-build.sh [--db $dbname]
	--db: directory to create database in (Default: $dbname)
	-j/--jar: location of GRIDSS jar
	-h/--help: display this message"
OPTIONS=hj:
LONGOPTS=help,db:,jar:
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
			dbname=$2
			shift 2
			;;
		-j|--jar)
			GRIDSS_JAR="$2"
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
if [[ "$@" != "" ]] ; then
	echo "$USAGE_MESSAGE" 1>&2
	exit $EX_USAGE
fi
write_status() {
	echo "$(date): $1" 1>&2
}
for tool in kraken2-build samtools gunzip wget awk tar dustmasker rsync java ; do
	if ! which $tool >/dev/null; then
		echo "Error: unable to find $tool on \$PATH" 1>&2
		exit $EX_CONFIG
	fi
	write_status "Found $(which $tool)"
done
### Find the jars
find_jar() {
	env_name=$1
	if [[ -f "${!env_name:-}" ]] ; then
		echo "${!env_name}"
	else
		write_status "Unable to find $2 jar. Specify using the environment variant $env_name, or the --jar command line parameter."
		exit $EX_NOINPUT
	fi
}
gridss_jar=$(find_jar GRIDSS_JAR gridss)
write_status "Using GRIDSS jar $gridss_jar"

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
for f in $(find $dbname/ -name '*.fna') ; do
	# TODO: CreateSequenceDictionary
	samtools faidx $f
	rm -f $f.dict
	java -cp $GRIDSS_JAR \
		picard.cmdline.PicardCommandLine \
		CreateSequenceDictionary \
		I=$f \
		R=$f.dict
done

cd $dbname
cd ..
tar -czvf virusbreakend.db.$(basename $dbname).tar.gz \
	$(basename $dbname)/*.k2d \
	$(basename $dbname)/taxonomy/nodes.dmp \
	$(basename $dbname)/virushostdb.tsv \
	$(basename $dbname)/library/viral/*.fna* \
	$(basename $dbname)/library/added/*.fna* \

write_status "VIRUSBreakend build successful"
write_status "The full build (including intermediate files) can be found in $dbname"
write_status "An archive containing only the files required by VIRUSBreakend can be found at $(realpath $dbname/../virusbreakend.db.$(basename $dbname).tar.gz)"

trap - EXIT
exit 0 # success!

