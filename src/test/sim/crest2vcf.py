#!/usr/local/bioinf/bin/python
#
# converts GASVPro output into pseudo-VCF
#
import sys

print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##ALT=<ID=DUP,Description="">
##ALT=<ID=DEL,Description="">
##ALT=<ID=INV,Description="">
##INFO=<ID=left_chr,Number=1,Type=String,Description="">
##INFO=<ID=left_pos,Number=1,Type=Integer,Description="">
##INFO=<ID=left_strand,Number=1,Type=String,Description="">
##INFO=<ID=left_softclipped_read_count,Number=1,Type=Integer,Description="">
##INFO=<ID=right_chr,Number=1,Type=String,Description="">
##INFO=<ID=right_pos,Number=1,Type=Integer,Description="">
##INFO=<ID=right_strand,Number=1,Type=String,Description="">
##INFO=<ID=right_softclipped_read_count,Number=1,Type=Integer,Description="">
##INFO=<ID=type,Number=1,Type=String,Description="">
##INFO=<ID=left_coverage,Number=1,Type=Integer,Description="">
##INFO=<ID=right_coverage,Number=1,Type=Integer,Description="">
##INFO=<ID=left_assembly_length,Number=1,Type=Integer,Description="assembled length at left_pos">
##INFO=<ID=right_assembly_length,Number=1,Type=Integer,Description="assembled length at right_pos">
##INFO=<ID=left_percent_identity,Number=1,Type=Float,Description="average percent identity at left_pos">
##INFO=<ID=left_percent_multimapping,Number=1,Type=Float,Description="percent of non-unique mapping reads at left_pos">
##INFO=<ID=right_percent_identity,Number=1,Type=Float,Description="average percent identity at right_pos">
##INFO=<ID=right_percent_multimapping,Number=1,Type=Float,Description="percent of non-unique mapping reads at right_pos">
##INFO=<ID=consensus_start_pos,Number=1,Type=Integer,Description="start position of consensus mapping to genome">
##INFO=<ID=consensus_start_chr,Number=1,Type=String,Description="starting chromosome of consensus mapping">
##INFO=<ID=consensus_mapping_start_pos,Number=1,Type=Integer,Description="position of the genomic mapping of consensus starting position">
##INFO=<ID=consensus_mapping_end_pos,Number=1,Type=Integer,Description="end position of consensus mapping to genome">
##INFO=<ID=consensus_end_chr,Number=1,Type=String,Description="ending chromsome of consnesus mapping">
##INFO=<ID=consensus_mapping_end_pos2,Number=1,Type=Integer,Description="position of genomic mapping of consensus ending posiiton">
##INFO=<ID=consensus,Number=1,Type=String,Description="consensus sequences">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""
"""The program will generate a *.predSV.txt file.  The filename will be the input
bam with .predSV.txt appended unless you specify the -p parameter.  Also the 
STDERR output has the full list of SVs, including rejected ones.  The output 
file *.predSV.txt has the following tab-delimited columns: left_chr, left_pos, 
left_strand, # of left soft-clipped reads, right_chr, right_pos, right_strand, 
# right soft-clipped reads, SV type, coverage at left_pos, coverage at 
right_pos, assembled length at left_pos, assembled length at right_pos,
average percent identity at left_pos, percent of non-unique mapping reads at
left_pos, average percent identity at right_pos, percent of non-unique mapping
reads at right_pos, start position of consensus mapping to genome,
starting chromosome of consensus mapping, position of the genomic mapping of
consensus starting position, end position of consensus mapping to genome,
ending chromsome of consnesus mapping, position of genomic mapping of
consensus ending posiiton, and consensus sequences.  For inversion(INV), the 
last 7 fields will be repeated to reflect the fact two different breakpoints 
are needed to identify an INV event.
"""

commonInfo=";left_chr={0};left_pos={1};left_strand={2};left_softclipped_read_count={3};right_chr={4};right_pos={5};right_strand={6};right_softclipped_read_count={7};left_coverage={9};right_coverage={10};left_assembly_length={11};right_assembly_length={12};left_percent_identity={13};left_percent_multimapping={14};right_percent_identity={15};right_percent_multimapping={16};consensus_start_pos={17};consensus_start_chr={18};consensus_mapping_start_pos={19};consensus_mapping_end_pos={20};consensus_end_chr={21};consensus_mapping_end_pos2={22};consensus={23}"

def toVcfBreakend(localChr, localPos, localPositive, remoteChr, remotePos, remotePositive):
	if remotePositive:
		remote = "]" + remoteChr + ":" + str(remotePos) + "]"
	else:
		remote = "[" + remoteChr + ":" + str(remotePos) + "["
	if localPositive:
		return "N" + remote
	else:
		return remote + "N"
	
i = 0
for line in sys.stdin:
	i += 1
	input = line.strip().split('	')
	left_chr = input[0]
	left_pos = int(input[1])
	left_strand = input[2]
	left_softclipped_read_count = int(input[3])
	right_chr = input[4]
	right_pos = int(input[5])
	right_strand = input[6]
	right_softclipped_read_count = int(input[7])
	svType = input[8]
	left_coverage = int(input[9])
	right_coverage = int(input[10])
	left_assembly_length = int(input[11]) #assembled length at left_pos
	right_assembly_length = int(input[12]) #assembled length at right_pos,
	left_percent_identity = float(input[13]) #average percent identity at left_pos
	left_percent_multimapping = float(input[14]) #percent of non-unique mapping reads at left_pos
	right_percent_identity = float(input[15]) # average percent identity at right_pos
	right_percent_multimapping = float(input[16]) # percent of non-unique mapping reads at right_pos
	consensus_start_pos = int(input[17]) # start position of consensus mapping to genome
	consensus_start_chr = input[18] # starting chromosome of consensus mapping
	consensus_mapping_start_pos = input[19] # position of the genomic mapping of consensus starting position
	consensus_mapping_end_pos = input[20] # end position of consensus mapping to genome
	consensus_end_chr = input[21] # ending chromsome of consnesus mapping
	consensus_mapping_end_pos = input[22] # position of genomic mapping of consensus ending posiiton (?same as input[20]?)
	consensus = input[23] # consensus sequences
	
	#chrX    55678965        +        13      chrX     52886736        +       0        INS     26      24       55      0    0.730555555555556    0.583333333333333       1        0       1       
	#chrX     55678880        141     chrX     52886790      ACATACTCTTTTGTCTTTGTCTTTATGCCCGTGTTCATCCTCCTTTGTTCAGTCCAGCAAGGTCTGCAGCATTATAAAGTTCAAAGGCATGGGAACCTAGAGCTGCCCCTTCTGTCTTTCTTTTAAGTAAGGTCCAAAGGT
	id = "line" + str(i)
	#http://www.stjuderesearch.org/site/lab/zhang
	if svType == "INS":
		# is a tandem duplication according to the CREST readme
		print("{0}	{5}	{24}	N	<DUP>	.	PASS	SVTYPE=DUP;SVLEN={25}" + commonInfo).format(
			left_chr, left_pos, left_strand, left_softclipped_read_count, # 0-3
			right_chr, right_pos, right_strand, right_softclipped_read_count, # 4-7
			svType, left_coverage, right_coverage, # 8-10
			left_assembly_length, right_assembly_length,
			left_percent_identity,left_percent_multimapping, right_percent_identity, right_percent_multimapping,
			consensus_start_pos, consensus_start_chr, consensus_mapping_start_pos, consensus_mapping_end_pos, consensus_end_chr, consensus_mapping_end_pos,
			consensus,
			id, #24
			left_pos - right_pos, #25
			)
	elif svType == "DEL":
		print("{0}	{1}	{24}	N	<DEL>	.	PASS	SVTYPE={8};SVLEN={25}" + commonInfo).format(
			left_chr, left_pos, left_strand, left_softclipped_read_count, # 0-3
			right_chr, right_pos, right_strand, right_softclipped_read_count, # 4-7
			svType, left_coverage, right_coverage, # 8-10
			left_assembly_length, right_assembly_length,
			left_percent_identity,left_percent_multimapping, right_percent_identity, right_percent_multimapping,
			consensus_start_pos, consensus_start_chr, consensus_mapping_start_pos, consensus_mapping_end_pos, consensus_end_chr, consensus_mapping_end_pos,
			consensus,
			id, #24
			right_pos - left_pos, #25
			)
	elif svType in ("CTX", "ITX"):
		print("{0}	{1}	{24}o	N	{25}	.	PASS	SVTYPE=BND;PARID={24}h;EVENT={8}{24}\n" +
			"{4}	{5}	{24}h	N	{26}	.	PASS	SVTYPE=BND;PARID={24}o;EVENT={8}{24}" + commonInfo).format(
			left_chr, left_pos, left_strand, left_softclipped_read_count, # 0-3
			right_chr, right_pos, right_strand, right_softclipped_read_count, # 4-7
			svType, left_coverage, right_coverage, # 8-10
			left_assembly_length, right_assembly_length, # 11-12
			left_percent_identity,left_percent_multimapping, right_percent_identity, right_percent_multimapping, # 13-16
			consensus_start_pos, consensus_start_chr, consensus_mapping_start_pos, consensus_mapping_end_pos, consensus_end_chr, consensus_mapping_end_pos, # 17-22
			consensus, #23
			id, # 24
			toVcfBreakend(left_chr, left_pos, left_strand=="+", right_chr, right_pos, right_strand!="+"), #25
			toVcfBreakend(right_chr, right_pos, right_strand!="+", left_chr, left_pos, left_strand=="+") #26
			)
	elif svType == "INV":
		# chr12   198894  +       20      chr12   202990  -       32      ITX     87      99      47      49      0.898709677419355       0       0.8875  0       1       chr12   202894  194     chr12   198798  AGAGTTGGAGGTGGGGCCTAACTCGAGGTATCTGGGTTATGGGAGTGGCTTCCTCACGAATAGATTAACACCCTCTCTGGGGAGTTGGGGGTGAGTGTAAATAAAGGAAGGTGTCCACAAAAGACTCCATATTGTATGATTCCACTTATATAATATTCTAGAAAAGACACAACTATATTGACAGAGGTAAGTTT
		print("{0}	{1}	{24}	N	<INV>	.	PASS	SVTYPE={8};SVLEN={25}" + commonInfo).format(
			left_chr, left_pos, left_strand, left_softclipped_read_count, # 0-3
			right_chr, right_pos, right_strand, right_softclipped_read_count, # 4-7
			svType, left_coverage, right_coverage, # 8-10
			left_assembly_length, right_assembly_length,
			left_percent_identity,left_percent_multimapping, right_percent_identity, right_percent_multimapping,
			consensus_start_pos, consensus_start_chr, consensus_mapping_start_pos, consensus_mapping_end_pos, consensus_end_chr, consensus_mapping_end_pos,
			consensus,
			id, #24
			right_pos - left_pos, #25
			)
	else:
		print "Undocumented variant type: " + line
	


