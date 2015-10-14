#!/usr/local/bioinf/bin/python
#
# converts breakdancer output into pseudo-VCF
#
import sys

fragSize = int(sys.argv[1])

print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##ALT=<ID=DEL,Description="">
##ALT=<ID=INS,Description="">
##ALT=<ID=INV,Description="">
##INFO=<ID=Chr1,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Pos1,Number=1,Type=Integer,Description="BreakDancer output">
##INFO=<ID=Orient1,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Chr2,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Pos2,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Orient2,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Type,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Size,Number=1,Type=Integer,Description="BreakDancer output">
##INFO=<ID=Score,Number=1,Type=Float,Description="BreakDancer output">
##INFO=<ID=num_Reads,Number=1,Type=Integer,Description="BreakDancer output">
##INFO=<ID=UNKNOWN_ORIENTATION,Number=0,Type=Flag,Description="Orientation of structural variation is unknown">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""
# CIPOS removed from output as rightPos appears to mean different things for different INV rows
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">


commonInfo=";Chr1={0};Pos1={1};Orient1={2};Chr2={3};Pos2={4};Orient2={5};Type={6};Size={7};Score={8};num_Reads={9}"

def toVcfBreakend(localChr, localPos, localPositive, remoteChr, remotePos, remotePositive):
	if remotePositive:
		remote = "]" + remoteChr + ":" + str(remotePos) + "]"
	else:
		remote = "[" + remoteChr + ":" + str(remotePos) + "["
	if localPositive:
		return "N" + remote
	else:
		return remote + "N"

def isRefToBreakend(orientation):
	str = orientation.replace("-", "+").split("+")
	plusCount = int(str[0])
	minusCount = int(str[1])
	return plusCount >= minusCount

i = 0
for line in sys.stdin:
	i = i + 1
	if line[0] == '#': continue
	#Chr1	Pos1	Orientation1	Chr2	Pos2	Orientation2	Type	Size	Score	num_Reads	num_Reads_lib	c1c173e1dd7a0e58ff9bd127c0cdce16.sc.bam
	id = "line" + str(i)
	input = line.split('	')
	leftChr = input[0]
	leftPos = int(input[1])
	leftOrientation = input[2]
	rightChr = input[3]
	rightPos = int(input[4])
	rightOrientation = input[5]
	type = input[6]
	size = int(input[7])
	score = int(input[8])
	numreads = int(input[9])
	if type == "ITX":
		print (
			"{0}	{1}	{10}a	N	N[{3}:{4}[	{8}	PASS	IMPRECISE;UNKNOWN_ORIENTATION;SVLEN={7};SVTYPE=BND;PARID={10}b;EVENT={6}{10}\n" +
			"{3}	{4}	{10}b	N	]{0}:{11}]N	{8}	PASS	IMPRECISE;UNKNOWN_ORIENTATION;SVLEN={7};SVTYPE=BND;PARID={10}a;EVENT={6}{10}\n" + 
			"{3}	{4}	{10}c	N	N[{0}:{1}[	{8}	PASS	IMPRECISE;UNKNOWN_ORIENTATION;SVLEN={7};SVTYPE=BND;PARID={10}d;EVENT={6}{10}\n" +
			"{0}	{11}	{10}d	N	]{3}:{4}]N	{8}	PASS	IMPRECISE;UNKNOWN_ORIENTATION;SVLEN={7};SVTYPE=BND;PARID={10}c;EVENT={6}{10}" +
			commonInfo
			).format(
			leftChr, leftPos, leftOrientation, rightChr, rightPos, rightOrientation, type, size, score, numreads, # 0-9
			id, #10
			leftPos + size) # 11
	elif type == "CTX":
		print (
			"{0}	{1}	{10}o	N	{11}	{8}	PASS	IMPRECISE;UNKNOWN_ORIENTATION;SVLEN={7};SVTYPE=BND;PARID={10}h;EVENT={6}{10}\n" +
			"{3}	{4}	{10}h	N	{12}	{8}	PASS	IMPRECISE;UNKNOWN_ORIENTATION;SVLEN={7};SVTYPE=BND;PARID={10}o;EVENT={6}{10}" +
			commonInfo).format(
			leftChr, leftPos, leftOrientation, rightChr, rightPos, rightOrientation, type, size, score, numreads, # 0-9
			id, #10
			toVcfBreakend(leftChr, leftPos, isRefToBreakend(leftOrientation), rightChr, rightPos, isRefToBreakend(rightOrientation)), #11
			toVcfBreakend(rightChr, rightPos, isRefToBreakend(rightOrientation), leftChr, leftPos, isRefToBreakend(leftOrientation))) # 12
	elif type == "INS":
		print (
			"{0}	{1}	{10}	N	<{6}>	{8}	PASS	IMPRECISE;END={1};SVLEN={11};SVTYPE={6}" + commonInfo # ;CIPOS=0,{12}
			).format(
			leftChr, leftPos, leftOrientation, rightChr, rightPos, rightOrientation, type, size, score, numreads, # 0-9
			id, #10
			-size, #11
			rightPos - leftPos #12 - error margin of call
			)
	elif type == "DEL":
		print (
			"{0}	{1}	{10}	N	<{6}>	{8}	PASS	IMPRECISE;END={14};SVLEN={11};SVTYPE={6}" + commonInfo # ;CIPOS={12},{13}
			).format(
			leftChr, leftPos, leftOrientation, rightChr, rightPos, rightOrientation, type, size, score, numreads, # 0-9
			id, #10
			-size, #11
			min(0, rightPos - leftPos - size), max(0, rightPos - leftPos - size), #12-13 - implied error margin of call
			leftPos + size #13 - end of variant at called position
			) 
	elif type == "INV":
		# what does rightPos mean for INV? these output lines don't make any sense at all
		#chr12   512749  0+2-    chr12   512812  0+2-    INV     74454768        41      2       0.sc.bam|2     NA
		#chr12   311335  41+0-   chr12   912971  41+0-   INV     601313  99      41      0.sc.bam|41    77.23
		
		# breakdancer systematically reports incorrect size and positions:
		# positions incorrect by 1 fragment length
		# length reported 1 fragment size too short
		leftPos -= fragSize
		rightPos -= fragSize
		size += fragSize
		print (
			#chr12	512749	0+2-	chr12	512812	0+2-	INV	74454768	41	2	/home/users/allstaff/cameron.d/i/data.fastcompare/3b57e681753bebdbde7eec2eb9da6510.sc.bam|2	NA
			"{0}	{1}	{10}	N	<{6}>	{8}	PASS	IMPRECISE;END={4};SVLEN={7};SVTYPE={6}" + commonInfo # ;CIPOS={12},{13}
			).format(
			leftChr, leftPos, leftOrientation, rightChr, rightPos, rightOrientation, type, size, score, numreads, # 0-9
			id, #10
			size #11
			)
	else:
		print "UNHANDLED VARIANT CALL" + line
	



