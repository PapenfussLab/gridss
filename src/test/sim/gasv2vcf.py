#!/usr/local/bioinf/bin/python
#
# converts GASVPro output into VCF
#
import sys

print """##fileformat=VCFv4.1
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around the length of the inserted material between breakends">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=LOCALISATION,Number=1,Type=Float,Description="Localisation">
##INFO=<ID=PRS,Number=1,Type=Float,Description="Num PRS">
##INFO=<ID=LLR,Number=1,Type=Float,Description="LogLikelihoodRatio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

lookup = {}
for line in open(sys.argv[1]):
	chrIndex = str.strip(line).split('	')
	lookup[chrIndex[1]] = chrIndex[0]

for line in sys.stdin:
	if line[0] == '#': continue
	#Cluster_ID:    LeftChr:        LeftBreakPoint: RightChr:       RightBreakPoint:        Num PRS:        Localization:   Type:   LogLikelihoodRatio:
	input = line.split('	')
	clusterid = input[0]
	leftChr = input[1]
	leftStart = int(input[2].split(',')[0])
	leftEnd = int(input[2].split(',')[1])
	rightChr = input[3]
	rightStart = int(input[4].split(',')[0])
	rightEnd = int(input[4].split(',')[1])
	numprs = input[5]
	localisation = float(input[6])
	type = input[7] 
	llr = float(input[8])
	
	startMidPos = (leftEnd + leftStart) / 2
	endMidPos = (rightEnd + rightStart) / 2
	svLen = endMidPos - startMidPos
	
	if type[0] == 'D':
		svType = "DEL"
		svLen = -1 * svLen
	if type[0] == 'I':
		svType = "INS"
	if leftChr == rightChr:
		print "{0}	{1}	{2}	.	.	{3}	PASS	IMPRECISE;END={6};SVLEN={5};SVTYPE={4};CIPOS={7},{8};CIEND={9},{10};LOCALISATION={11};PRS={12};LLR={13}".format(
			lookup[leftChr],
			startMidPos,
			clusterid,
			numprs,
			svType,
			svLen,
			endMidPos,
			leftStart - startMidPos,
			leftEnd - startMidPos,
			rightStart - endMidPos,
			rightEnd - endMidPos,
			localisation,
			numprs,
			llr
			)
	else:
		print "Unhandled interchromsomal event: " + input
		
# Output format before convertClusters has been run
# for line in sys.stdin:
	# if line[0] == '#': continue
	# #Cluster_ID:	Num PRS:	Localization:	Type:	List of PRS:	LeftChr:	RightChr:	Boundary Points:	Log_Likelihood_Ratio:
	# input = line.split('	')
	# clusterid = input[0]
	# pairCount = input[1]
	# localisation = float(input[2])
	# type = input[3]
	# pairList = input[4]
	# leftChr = input[5]
	# rightChr = input[6]
	# boundaryPoints = map(int, input[7].split(','))
	# logLikelihood = float(input[8])
	# # transform
	# points = sorted(boundaryPoints)
	# lenp = len(points)
	# halflenp = lenp / 2
	# startMidPos = sum(points[0:halflenp]) / halflenp
	# endMidPos = sum(points[halflenp:lenp]) / halflenp
	# svLen = endMidPos - startMidPos
	# if type[0] == 'D':
		# svType = "DEL"
		# svLen = -1 * svLen
	# if type[0] == 'I':
		# svType = "INS"
	# if leftChr == rightChr:
		# print "{0}	{1}	{2}	.	.	{3}	PASS	IMPRECISE;END={6};SVLEN={5};SVTYPE={4};CIPOS={7},{8};CIEND={9},{10};CILEN={11},{12};LOCALISATION={13}".format(
			# lookup[leftChr],
			# startMidPos,
			# clusterid,
			# logLikelihood / 10,
			# svType,
			# svLen,
			# endMidPos,
			# points[0] - startMidPos,
			# points[halflenp-1] - startMidPos,
			# points[halflenp] - endMidPos,
			# points[lenp-1] - endMidPos,
			# points[halflenp] - points[halflenp-1],
			# points[lenp-1] - points[0],
			# localisation
			# )
	# else:
		# print "Unhandled interchromsomal event: " + input
		

