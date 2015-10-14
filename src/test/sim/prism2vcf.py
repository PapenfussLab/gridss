#!/usr/local/bioinf/bin/python
#
# converts breakdancer output into pseudo-VCF
#
import sys

print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

i = 0
for line in sys.stdin:
	i = i + 1
	if line[0] == '@': continue
	#@chr	gap_start	gap_end	gap_size	type	support_read_num	contig_match	contig_snp	ab_isize_0	ab_isize_1	avg_isize	mismatch_1	mismatch_2	mismatch_3	left_match	right_match	score1	score2	score3	mate_soft_clipped	mate_unmapped	sv_seq
	input = line.split('	')
	chr = input[0]
	gap_start = input[1]
	gap_end = input[2]
	gap_size = input[3]
	type = input[4]
	support_read_num = input[5]
	contig_match = input[6]
	contig_snp = input[7]
	ab_isize_0 = input[8]
	ab_isize_1 = input[9]
	avg_isize = input[10]
	mismatch_1 = input[11]
	mismatch_2 = input[12]
	mismatch_3 = input[13]
	left_match = input[14]
	right_match = input[15]
	score1 = input[16]
	score2 = input[17]
	score3 = input[18]
	mate_soft_clipped = input[19]
	mate_unmapped = input[20]
	sv_seq = input[21]
	
	size = int(gap_size)
	if type == "DEL":
		size = size * -1
	print "{0}	{1}	{2}	N	<{6}>	{3}	PASS	END={4};SVLEN={5};SVTYPE={6}".format(
		chr,
		gap_start,
		".",
		support_read_num,
		gap_end,
		size,
		type
		)

