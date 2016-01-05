#!/usr/local/bioinf/bin/python
# Daniel Cameron
#
# converts socrates paired output to VCF breakpoints
#
import sys
import os

print """##fileformat=VCFv4.1
##FILTER=<ID=UNPAIRED,Description="One one breakend has soft clip support">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=ANCHCONS,Number=1,Type=String,Description="Anchor consensus sequence">
##INFO=<ID=REALNCONS,Number=1,Type=String,Description="Consensus sequence of realigned soft-clipped reads aligned to mate end">
##INFO=<ID=NLSC,Number=1,Type=Integer,Description="Number of long soft-clips supporting this breakpoint">
##INFO=<ID=NSSC,Number=1,Type=Integer,Description="Number of short soft-clips supporting this breakpoint">
##INFO=<ID=BLSC,Number=1,Type=Integer,Description="Total number of long soft-clips bases supporting this breakpoint">
##INFO=<ID=BSSC,Number=1,Type=Integer,Description="Total number of short soft-clips bases supporting this breakpoint">
##INFO=<ID=LSSC,Number=1,Type=Integer,Description="Length of longest sort soft-clip supporting this breakpoint">
##INFO=<ID=COMMENT,Number=1,Type=Integer,Description="Undefined comment field">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

#C1_realign	C1_realign_dir	C1_realign_consensus	C1_anchor	C1_anchor_dir	C1_anchor_consensus	C1_long_support	C1_long_support_bases	C1_short_support	C1_short_support_bases	C1_short_support_max_len	C1_avg_realign_mapq	C2_realign	C2_realign_dir	C2_realign_consensus	C2_anchor	C2_anchor_dir	C2_anchor_consensus	C2_long_support	C2_long_support_bases	C2_short_support	C2_short_support_bases	C2_short_support_max_len	C2_avg_realign_mapq	BP_condition
'''
There are two files generated: the paired and unpiared outputs. The paired output contains the best results, while the unpaired contains any soft-clips that have realigned anywhere else in the process re-alignments stage of the algorithm. The unpaired results are included for completeness as they can contain some useful information, but overall this output is comprised of false positives due to mapping errors and other artefacts.

The paired output contains various columns of information that describe the location of the break point and the level of support:
	1. C1_realign - this genomic locus describes the position of the realigned soft-clips of cluster The position is a consensus if the realignements are not exactly the same for all soft-clips. The locus is that of the first soft- clipped base, so the one immediately next to the breakpoint.
	2. C1_realign_dir - this field takes either "+" or "-" and indicates whether the realigned softclips map upstream of C1_realign (+) or downstream (-) -- with respect to the reference genome.
	3. C1_realign_consensus - a consensus sequence made of soft-clips in cluster 1. An asterisk optionally marks the position of the breakpoint within the consensus (if the position was unanimous).
	4. C1_anchor - a second locus describing the anchor region of the cluster. The anchor is defined by the reads that were mapped in the initial alignments and which soft-clips formed the earlier columns. The postion is the consensus of positions before the first soft-clipped base.
	5. C1_anchor_dir - analogous to above this field describes whether the anchor region is upstream ("+") or downstream ("-") of the breakpoint in the reference.
	6. C1_anchor_consensus - the consensus sequence of the anchor reads.
	7. C1_long_support - this number counts the number of "long" soft-clips (as specified during the run of Socrates) that support the cluster C1. This number is at least 1, as there would not be a cluster without a realigned soft-clip.
	8. C1_long_support_bases - the number of nucleotides in the long support of the preceding column counted and reported in this column.
	9. C1_short_support - similarly, the short support is counted in number...
	10. C1_short_support_bases - ... and nucleotides.
	11. C1_avg_realign_mapq - this last column for C1 summarizes the average mapping quality of the anchor reads. It can be a helpful filter criteria and a minimum value can be specified at the launch of Socrates.
	12-22. C2 columns - all columns described above are repeated for C2. There is only one noticeable difference: C2 can be a cluster formed without realigned soft-clips (see "short SC cluster" below), which leads to empty consensus sequences (there are double-tabs in the output, which can be quite nasty to deal with. Apologies!) and long-support values set to 0. Finally, the column labelled "BP_condition" describes the nature of the fusion event. It can take any of five values:
		1. Blunt-end joining: the most straight forward case of a clean join (none of the below).
		2. Micro-homology: Xbp homology found! (XXX): the two joined regions are identical for X bases across the break. Therefore the true location of the breakpoint is only known within those boundaries.
		3. Inserted sequence: XXX: There is a short bit of sequence inserted in between the two loci of the fusion. The sequence is either untemplated or from somewhere else in the genome (but too short to map).
		4. unequal distances of realigned breakpoint to anchor breakpoint: X v Y: In this case the realignment and anchor loci of the two paired clusters do not support the exact same coordinate for a fusion (|X-Y| indicates the difference). This is usually due to mis-mappings.
		5. Unequal inserted sequence: XXX v Y: an insert occurs as above, but Socrates was unable to determine the exact sequence. One of the values should contain the correct sequence.
		6. short SC cluster: In most experiments the most prevalent type, yet the least trustworthy. Only one side of the breakpoint is supported by realiged split reads, the other by short bits of soft-clipped sequence only. This sort of cluster pairing makes Socrates very sensitive, but introduces false positives.
'''

def toVcfBreakend(localChr, localPos, localDir, remoteChr, remotePos, remoteDir, untemplated):
	if remoteDir == "+":
		remote = "]" + remoteChr + ":" + str(remotePos) + "]"
	else:
		remote = "[" + remoteChr + ":" + str(remotePos) + "["
	if localDir == "+":
		return "N" + untemplated + remote
	else:
		return remote + untemplated + "N"

def toVcfPos(pos, dir):
	if dir == "-":
		return pos + 1
	else:
		return pos - 1
		
comp_lookup = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
			  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n',} 
			  
def revComp(s):
	return "".join([comp_lookup[base] for base in list(s)])

def processPairedLine(i, line):
	input = line.split('	')
	C1_realign = input[0]
	C1_realign_dir = input[1]
	C1_realign_consensus = input[2]
	C1_anchor = input[3]
	C1_anchor_dir = input[4]
	C1_anchor_consensus = input[5]
	C1_long_support = int(input[6])
	C1_long_support_bases = int(input[7])
	C1_short_support = int(input[8])
	C1_short_support_bases = int(input[9])
	C1_short_support_max_len = int(input[10])
	C1_avg_realign_mapq = float(input[11])
	C2_realign = input[12]
	C2_realign_dir = input[13]
	C2_realign_consensus= input[14]
	C2_anchor = input[15]
	C2_anchor_dir = input[16]
	C2_anchor_consensus = input[17]
	C2_long_support = int(input[18])
	C2_long_support_bases = int(input[19])
	C2_short_support = int(input[20])
	C2_short_support_bases = int(input[21])
	C2_short_support_max_len = int(input[22])
	C2_avg_realign_mapq = float(input[23])
	comment = input[24]
	comment = comment.strip()
	comment = comment.replace(" ", "").strip()
	
	# calculate
	C1_realign_chr = C1_realign.split(":")[0]
	C1_realign_pos = int(C1_realign.split(":")[1])
	C2_realign_chr = C2_realign.split(":")[0]
	C2_realign_pos = int(C2_realign.split(":")[1])
	C1_anchor_chr = C1_anchor.split(":")[0]
	C1_anchor_pos = int(C1_anchor.split(":")[1])
	C2_anchor_chr = C2_anchor.split(":")[0]
	C2_anchor_pos = int(C2_anchor.split(":")[1])
	min_avg_mapq = min(C1_avg_realign_mapq, C2_avg_realign_mapq)
	min_long_support = min(C1_long_support, C2_long_support)
	
	# quality score proxy
	score = C1_long_support + C1_short_support + C2_long_support + C2_short_support
	
	# sanity check
	if C1_realign_chr != C2_anchor_chr:
		print >> sys.stderr, "Sanity check failure on line {0}: C1_realign_chr of {1} does not match C2_anchor_chr of {2}".format(
			i, C1_realign_chr, C2_anchor_chr)
		return
	if C2_realign_chr != C1_anchor_chr:
		print >> sys.stderr, "Sanity check failure on line {0}: C2_realign_chr of {1} does not match C1_anchor_chr of {2}".format(
			i, C2_realign_chr, C1_anchor_chr)
		return
	if C1_realign_dir != C2_anchor_dir:
		print >> sys.stderr, "Sanity check failure on line {0}: C1_realign_dir of {1} does not match C2_anchor_dir of {2}".format(
			i, C1_realign_dir, C2_anchor_dir)
		return
	if C2_realign_dir != C1_anchor_dir:
		print >> sys.stderr, "Sanity check failure on line {0}: C2_realign_dir of {1} does not match C1_anchor_dir of {2}".format(
			i, C2_realign_dir, C1_anchor_dir)
		return
	
	# calculate mapping offsets to determine microhomology or untemplated sequence
	C1_microhomology_len = C2_realign_pos - C1_anchor_pos
	if C2_realign_dir == "+" : 
		C1_microhomology_len *= -1
	C2_microhomology_len = C1_realign_pos - C2_anchor_pos
	if C1_realign_dir == "+" : 
		C2_microhomology_len *= -1
	
	# unequal distances of realigned breakpoint to anchor breakpoint
	# we'll just use C1 for the VCF breakpoint call
#	if C1_microhomology_len != C2_microhomology_len:
#		print >> sys.stderr, "Sanity check failure on line {0}: C1_microhomology_len of {1} does not match C2_microhomology_len of {2}".format(
#			i, C1_microhomology_len, C2_microhomology_len)
#		continue

	C1_cipos = ""
	C2_cipos = ""
	if C1_microhomology_len > 0:
		# C1 positioning is used as the arbitrary breakpoint
		# microhomology
		C1_cipos = "CIPOS=-{0},0;".format(C1_microhomology_len)
		C2_cipos = "CIPOS=0,{0};".format(C1_microhomology_len)
	
	C1_untemplatedSequence = ""
	C1_realign_bases = C1_realign_consensus
	if "*" in C1_realign_consensus:
		C1_untemplatedSequence = C1_realign_consensus[:C1_realign_consensus.find("*")]
		C1_realign_bases = C1_realign_consensus[(C1_realign_consensus.find("*")+1):]
	
	C2_untemplatedSequence = C1_untemplatedSequence
	if C1_realign_dir == C2_anchor_dir:
		C2_untemplatedSequence = revComp(C1_untemplatedSequence)
	
	# adjust positions to first anchored based, instead of first breakend base
	C1_realign_pos = toVcfPos(C1_realign_pos, C1_realign_dir)
	C2_realign_pos = toVcfPos(C2_realign_pos, C2_realign_dir)
	C1_anchor_pos = toVcfPos(C1_anchor_pos, C1_anchor_dir)
	C2_anchor_pos = toVcfPos(C2_anchor_pos, C2_anchor_dir)
	
	C1_breakpoint = toVcfBreakend(C1_realign_chr, C1_realign_pos, C1_realign_dir, C1_anchor_chr, C1_anchor_pos, C1_anchor_dir, C1_untemplatedSequence)
	C2_breakpoint = toVcfBreakend(C1_anchor_chr, C1_anchor_pos, C1_anchor_dir, C1_realign_chr, C1_realign_pos, C1_realign_dir, C2_untemplatedSequence)
	
	# TODO offset call position to first anchored base
	# TODO work out what +- means for Socrates - opposite of what I expected!
	print "{1}	{2}	paired{0}o	N	{3}	{4}	PASS	{5}SVTYPE=BND;MATEID=paired{0}h;ANCHCONS={6};REALNCONS={7};NLSC={8};NSSC={9};BLSC={10};BSSC={11};LSSC={12};COMMENT={13}".format(
		i, C1_realign_chr, C1_realign_pos, C1_breakpoint,
		score,
		C1_cipos,
		C1_anchor_consensus, C2_realign_consensus,
		C1_long_support, C1_long_support_bases, C1_short_support, C1_short_support_bases, C1_short_support_max_len,
		comment.strip()
		)
	print "{1}	{2}	paired{0}h	N	{3}	{4}	PASS	{5}SVTYPE=BND;MATEID=paired{0}o;ANCHCONS={6};REALNCONS={7};NLSC={8};NSSC={9};BLSC={10};BSSC={11};LSSC={12};COMMENT={13}".format(
		i, C1_anchor_chr, C1_anchor_pos, C2_breakpoint,
		score,
		C2_cipos,
		C2_anchor_consensus, C2_realign_consensus,
		C2_long_support, C2_long_support_bases, C2_short_support, C2_short_support_bases, C2_short_support_max_len,
		comment.strip()
		)

def processUnpairedLine(i, line):
	input = line.split('	')
	C1_realign = input[0]
	C1_realign_dir = input[1]
	C1_anchor = input[3]
	C1_anchor_dir = input[4]
	C1_long_support = int(input[6])
	C1_long_support_bases = int(input[7])
	C1_short_support = int(input[8])
	C1_short_support_bases = int(input[9])
	C1_short_support_max_len = int(input[10])
	C1_avg_realign_mapq = float(input[11])
	
	# calculate
	C1_realign_chr = C1_realign.split(":")[0]
	C1_realign_pos = int(C1_realign.split(":")[1])
	C1_anchor_chr = C1_anchor.split(":")[0]
	C1_anchor_pos = int(C1_anchor.split(":")[1])
	
	# quality score proxy
	score = C1_long_support + C1_short_support
	

	# adjust positions to first anchored based, instead of first breakend base
	C1_realign_pos = toVcfPos(C1_realign_pos, C1_realign_dir)
	C1_anchor_pos = toVcfPos(C1_anchor_pos, C1_anchor_dir)
	
	C1_breakpoint = toVcfBreakend(C1_realign_chr, C1_realign_pos, C1_realign_dir, C1_anchor_chr, C1_anchor_pos, C1_anchor_dir, "")
	C2_breakpoint = toVcfBreakend(C1_anchor_chr, C1_anchor_pos, C1_anchor_dir, C1_realign_chr, C1_realign_pos, C1_realign_dir, "")
	
	print "{1}	{2}	unpaired{0}o	N	{3}	{4}	UNPAIRED	SVTYPE=BND;MATEID=unpaired{0}h".format(
		i, C1_realign_chr, C1_realign_pos, C1_breakpoint,
		score,
		)
	print "{1}	{2}	unpaired{0}h	N	{3}	{4}	UNPAIRED	SVTYPE=BND;MATEID=unpaired{0}o;NLSC={5};NSSC={6};BLSC={7};BSSC={8};LSSC={9}".format(
		i, C1_anchor_chr, C1_anchor_pos, C2_breakpoint,
		score,
		C1_long_support, C1_long_support_bases, C1_short_support, C1_short_support_bases, C1_short_support_max_len
		)


i = 0
for line in open( sys.argv[1] ):
	i += 1
	if line[0] == '#':
		# skip header line
		continue
	processPairedLine(i, line)

for line in open( sys.argv[2] ):
	i += 1
	if line[0] == '#':
		# skip header line
		continue
	processUnpairedLine(i, line)
