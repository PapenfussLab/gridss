conda activate kraken2
export PATH=~/dev/kraken2/kraken2-2.0.9-beta/bin:$PATH

# dbname parameter
# check for: kraken2-build samtools gunzip wget awk tar 
dbname=virusbreakenddb
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

