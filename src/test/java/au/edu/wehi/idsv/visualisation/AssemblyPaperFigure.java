package au.edu.wehi.idsv.visualisation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakpointFastqEncoding;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipEvidence;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

public class AssemblyPaperFigure extends IntermediateFilesTest {
	@Test
	public void generateFigureData() throws IOException {
		BreakpointSummary bs = new BreakpointSummary(2, FWD, 100, 100, 2, BWD, 5000, 5000);
		int k = 13;
		int readLength = 25;
		int fragSize = 75;
		ProcessingContext pc = getCommandlineContext(false);
		
		pc.getAssemblyParameters().k = k;
		pc.getConfig().getVisualisation().directory = new File(super.testFolder.getRoot().getAbsoluteFile(), "visualisation");
		pc.getConfig().getVisualisation().assembly = true;
		pc.getConfig().getVisualisation().assemblyProgress = true;
		pc.getConfig().getVisualisation().directory.mkdir();
				
		// Default 
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		List<SAMRecord> splitReads = new ArrayList<SAMRecord>();
		// lots of reads to set the expected fragment size
		for (int i = 1; i <= 1000; i++) {
			reads.add(RP(0, i, i + fragSize - readLength, readLength)[0]);
			reads.add(RP(0, i, i + fragSize - readLength, readLength)[1]);
		}
		// breakpoint from random 100 to random 
		// DPs
		addDP(reads, bs, readLength,  0, fragSize);
		addDP(reads, bs, readLength,  5, fragSize + 2);
		addDP(reads, bs, readLength, 20, fragSize - 4);
		// OEA
		addOEA(reads, bs, readLength, fragSize - readLength - 5, fragSize);
		addOEA(reads, bs, readLength, fragSize - readLength - 5, fragSize);
		
		addSR(reads, splitReads, bs, readLength, 10, 0);
		addSR(reads, splitReads, bs, readLength, 5, 0);
		addSR(reads, splitReads, bs, readLength, 13, 0);
		addSR(reads, splitReads, bs, readLength, 15, 0);
		
		createInput(reads);
		Collections.sort(splitReads, new Ordering<SAMRecord>() {
			@Override
			public int compare(SAMRecord left, SAMRecord right) {
				return ComparisonChain.start()
				.compare(BreakpointFastqEncoding.getEncodedReferenceIndex(left.getReadName()), BreakpointFastqEncoding.getEncodedReferenceIndex(right.getReadName()))
				.compare(BreakpointFastqEncoding.getEncodedStartPosition(left.getReadName()), BreakpointFastqEncoding.getEncodedStartPosition(right.getReadName()))
				.result();
			}
		});
		// assemble evidence
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(pc.getFileSystemContext().getRealignmentBam(input, 0), SortOrder.unsorted, splitReads.toArray(new SAMRecord[0]));
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), output);
		aes.ensureAssembled();
		// done
	}		
	private void addSR(List<SAMRecord> reads, List<SAMRecord> splitReads, BreakpointSummary bs, int readLength, int splitLength, int localledOvermappedBases) {
		int localLength = readLength - splitLength;
		assert(localledOvermappedBases == 0);
		SAMRecord r = Read(bs.referenceIndex, getStartBasePos(bs, localLength), bs.direction == FWD ? String.format("%dM%dS", localLength, splitLength) : String.format("%dS%dM", splitLength, localLength));
		String localBases = S(getRef(bs.referenceIndex, localLength, getStartBasePos(bs, localLength), bs.direction, localLength));
		String remoteBases = S(getRef(bs.remoteBreakend().referenceIndex, splitLength, getStartBasePos(bs.remoteBreakend(), splitLength), bs.remoteBreakend().direction, splitLength));
		if (bs.direction == FWD) r.setReadBases(B(localBases + remoteBases)); 
		else r.setReadBases(B(remoteBases + localBases));
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(SoftClipEvidence.create(SES(), bs.direction, r, null));
		assert(splitLength == fq.length());
		SAMRecord realign = Read(bs.remoteBreakend().referenceIndex, getStartBasePos(bs.remoteBreakend(), splitLength), String.format("%dM", splitLength));
		realign.setReadName(fq.getReadHeader());
		reads.add(r);
		splitReads.add(realign);
	}
	/**
	 * Adds discordant read pairs to the evidence set
	 * @param reads output list
	 * @param readLength
	 * @param offset offset from bp, zero indicates that the read aligns right up to the breakpoint
	 * @param fragSize
	 */
	private void addDP(List<SAMRecord> reads, BreakpointSummary bs, int readLength, int offset, int fragSize) {
		assert(bs.referenceIndex == bs.referenceIndex2);
		int fragBasesAtLocal = readLength + offset;
		int fragBasesAtRemote = fragSize - fragBasesAtLocal;
		assert(fragBasesAtLocal >= readLength);
		assert(fragBasesAtRemote >= readLength);
		assert(bs.direction == FWD);
		SAMRecord[] dp = DP(bs.referenceIndex, getStartBasePos(bs, readLength + offset), String.format("%dM", readLength), bs.direction == FWD,
				bs.remoteBreakend().referenceIndex, getStartBasePos(bs.remoteBreakend(), fragBasesAtRemote - (readLength - 1)), String.format("%dM", readLength), bs.remoteBreakend().direction == FWD);
		dp[0].setReadBases(B(S(RANDOM).substring(dp[0].getAlignmentStart() - 1,  dp[0].getAlignmentEnd())));
		dp[1].setReadBases(B(S(RANDOM).substring(dp[1].getAlignmentStart() - 1,  dp[1].getAlignmentEnd())));
		dp[0].setProperPairFlag(false);
		dp[1].setProperPairFlag(false);
		
		reads.add(dp[0]);
		reads.add(dp[1]);
	}
	/**
	 * Adds discordant read pairs to the evidence set
	 * @param reads output list
	 * @param readLength
	 * @param offset offset from bp, zero indicates that the read aligns right up to the breakpoint
	 * @param fragSize
	 */
	private void addOEA(List<SAMRecord> reads, BreakpointSummary bs, int readLength, int offset, int fragSize) {
		assert(bs.referenceIndex == bs.referenceIndex2);
		int fragBasesAtLocal = readLength + offset;
		int fragBasesAtRemote = fragSize - fragBasesAtLocal;
		assert(fragBasesAtLocal >= readLength);
		assert(fragBasesAtRemote > 0 && fragBasesAtRemote < readLength);
		assert(bs.referenceIndex == 2);
		int localSplitBases = readLength - fragBasesAtRemote; 
		SAMRecord[] dp = OEA(bs.referenceIndex, getStartBasePos(bs, readLength + offset), String.format("%dM", readLength), bs.direction == FWD);
		dp[0].setReadBases(B(S(RANDOM).substring(dp[0].getAlignmentStart() - 1,  dp[0].getAlignmentEnd())));
		// TODO: reverseComp this?
		dp[1].setReadBases(B(S(getRef(bs.referenceIndex, localSplitBases, bs.start, bs.direction, localSplitBases))
				+ S(getRef(bs.remoteBreakend().referenceIndex, fragBasesAtRemote, bs.remoteBreakend().start, bs.remoteBreakend().direction, fragBasesAtRemote))));
		reads.add(dp[0]);
		reads.add(dp[1]);
	}
	//@Test
	public void generateMissassemblyGraph() throws IOException {
		// [1,10] GTAAAAAC
		//fail();
	}
}
