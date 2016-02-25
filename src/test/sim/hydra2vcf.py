#!/usr/local/bioinf/bin/python
#
# converts breakdancer output into pseudo-VCF
#
import sys

print """##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=numDistinctPairs,Number=1,Type=Integer,Description="Number of distinct pairs in breakpoint. ">
##INFO=<ID=meanEditDist1,Number=1,Type=Float,Description="Mean edit distance observed in end1 of the breakpoint pairs. ">
##INFO=<ID=meanEditDist2,Number=1,Type=Float,Description="Mean edit distance observed in end2 of the breakpoint pairs. ">
##INFO=<ID=meanMappings1,Number=1,Type=Integer,Description="Mean number of mappings for end1 of all pairs in the breakpoint. ">
##INFO=<ID=meanMappings2,Number=1,Type=Integer,Description="Mean number of mappings for end2 of all pairs in the breakpoint. ">
##INFO=<ID=breakpointSize,Number=1,Type=Integer,Description="breakpointSize Size of the breakpoint. ">
##INFO=<ID=numMappings,Number=1,Type=Integer,Description="Total number of mappings included in the breakpoint. ">
##INFO=<ID=allWeightedSupport,Number=1,Type=Float,Description="Amount of weighted support from the mappings in the breakpoint. ">
##INFO=<ID=finalSupport,Number=1,Type=Float,Description="Amount of final support from the mappings in the breakpoint. ">
##INFO=<ID=finalWeightedSupport,Number=1,Type=Float,Description="Amount of final weighted support from the mappings in the breakpoint. ">
##INFO=<ID=numUniquePairs,Number=1,Type=Integer,Description="Number of pairs in the breakpoint that were uniquely mapped to the genome. ">
##INFO=<ID=numAnchoredPairs,Number=1,Type=Float,Description="Number of pairs in the breakpoint that were mapped to the genome in an "anchored" fashion (i.e. 1xN). ">
##INFO=<ID=numMultiplyMappedPairs,Number=1,Type=Integer,Description="Number of pairs in the breakpoint that were multiply mapped to the genome in fashion (i.e. NxN).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

def toVcfBreakend(localChr, localPos, localDir, remoteChr, remotePos, remoteDir):
	if remoteDir == "+":
		remote = "]" + remoteChr + ":" + str(remotePos) + "]"
	else:
		remote = "[" + remoteChr + ":" + str(remotePos) + "["
	if localDir == "+":
		return "N" + remote
	else:
		return remote + "N"

def toVcfPos(pos, dir):
	if dir == "-":
		return pos + 1
	else:
		return pos - 1

i = 0
for line in sys.stdin:
	i = i + 1
	if line[0] == '#': continue
	# https://code.google.com/archive/p/hydra-sv/wikis/FileFormats.wiki
	id = "line" + str(i)
	input = line.split('	')
	chrom1 = input[0] # 1. Chromosome for end 1 of the breakpoint
	start1 = int(input[1]) # 2. start1 Start position for end 1 of the breakpoint. 
	end1 = int(input[2]) # 3. end1 End position for end 1 of the breakpoint. 
	chrom2 = input[3] # 4. chrom2 Chromosome for end 2 of the breakpoint. 
	start2 = int(input[4]) # 5. start2 Start position for end 2 of the breakpoint.
	end2 = int(input[5]) # 6. end2 End position for end 2 of the breakpoint.
	breakpointId = input[6] # 7. breakpointId Unique Hydra breakpoint identifier.
	numDistinctPairs = int(input[7]) # 8. numDistinctPairs Number of distinct pairs in breakpoint. 
	strand1 = input[8] # 9. strand1 Orientation for the first end of the breakpoint. 
	strand2 = input[9] # 10. strand2 Orientation for the second end of the breakpoint. 
	meanEditDist1 = float(input[10]) # 11. meanEditDist1 Mean edit distance observed in end1 of the breakpoint pairs. 
	meanEditDist2 = float(input[11]) # 12. meanEditDist2 Mean edit distance observed in end2 of the breakpoint pairs. 
	meanMappings1 = float(input[12]) # 13. meanMappings1 Mean number of mappings for end1 of all pairs in the breakpoint. 
	meanMappings2 = float(input[13]) # 14. meanMappings2 Mean number of mappings for end2 of all pairs in the breakpoint. 
	breakpointSize = int(input[14]) # 15. breakpointSize Size of the breakpoint. 
	numMappings = int(input[15])# 16. numMappings Total number of mappings included in the breakpoint. 
	allWeightedSupport = float(input[16]) # 17. allWeightedSupport Amount of weighted support from the mappings in the breakpoint. 
	finalSupport = float(input[17]) # 18. finalSupport Amount of final support from the mappings in the breakpoint. 
	finalWeightedSupport = float(input[18]) # 19. finalWeightedSupport Amount of final weighted support from the mappings in the breakpoint. 
	numUniquePairs = int(input[19]) # 20. numUniquePairs Number of pairs in the breakpoint that were uniquely mapped to the genome. 
	numAnchoredPairs = float(input[20]) # 21. numAnchoredPairs Number of pairs in the breakpoint that were mapped to the genome in an "anchored" fashion (i.e. 1xN). 
	numMultiplyMappedPairs = int(input[21]) # 22. numMultiplyMappedPairs Number of pairs in the breakpoint that were multiply mapped to the genome in fashion (i.e. NxN).
	
	print ("{0}	{1}	hydra{2}o	N	{3}	{4}	PASS	SVTYPE=BND" +
		";MATEID=hydra{2}h;numDistinctPairs={5};meanEditDist1={6};meanEditDist2={7};meanMappings1={8};meanMappings2={9}{10}" + 
		";numMappings={11};allWeightedSupport={12};finalSupport={13};finalWeightedSupport={14};numUniquePairs={15};numAnchoredPairs={16};numMultiplyMappedPairs={17}"
		).format(
		chrom1, start1, breakpointId, finalSupport,
		toVcfBreakend(chrom1, start1, strand1, chrom2, start2, strand2),
		numDistinctPairs, meanEditDist1, meanEditDist2, meanMappings1, meanMappings2,
		"",#;SVLEN= #breakpointSize, # TODO convert to SVLEN
		numMappings, allWeightedSupport, finalSupport, finalWeightedSupport, numUniquePairs, numAnchoredPairs, numMultiplyMappedPairs
		)



