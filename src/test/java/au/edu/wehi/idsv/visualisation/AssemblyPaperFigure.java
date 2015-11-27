package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakpointFastqEncoding;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;

public class AssemblyPaperFigure extends IntermediateFilesTest {
	//@Test
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
		addDP(reads, null, bs, readLength,  0, fragSize, 0);
		addDP(reads, null, bs, readLength,  5, fragSize + 2, 0);
		addDP(reads, null, bs, readLength, 20, fragSize - 4, 0);
		// OEA
		addOEA(reads, bs, readLength, fragSize - readLength - 5, fragSize);
		addOEA(reads, bs, readLength, fragSize - readLength - 5, fragSize);
		
		addSR(reads, splitReads, null, bs, readLength, 10, 0);
		addSR(reads, splitReads, null, bs, readLength, 5, 0);
		addSR(reads, splitReads, null, bs, readLength, 13, 0);
		addSR(reads, splitReads, null, bs, readLength, 15, 0);
		
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
	private void addSR(List<SAMRecord> reads, List<SAMRecord> splitReads, List<SAMRecord> fullSplitReads, BreakpointSummary bs, int readLength, int splitLength, int localledOvermappedBases) {
		int localLength = readLength - splitLength;
		assert(localledOvermappedBases == 0);
		SAMRecord r = sr(bs, localLength, splitLength); 
		FastqRecord fq = BreakpointFastqEncoding.getRealignmentFastq(SoftClipEvidence.create(SES(), bs.direction, r, null)); 
		assert(splitLength == fq.length());
		SAMRecord realign = Read(bs.remoteBreakend().referenceIndex, getLeftmostBasePos(bs.remoteBreakend(), splitLength, 0), String.format("%dM", splitLength));
		realign.setReadName(fq.getReadHeader());
		realign.setReadBases(getRef(bs.referenceIndex2, splitLength, bs.start2, bs.direction2));
		reads.add(r);
		splitReads.add(realign);
		if (fullSplitReads != null) {
			SAMRecord remote = sr(bs.remoteBreakpoint(), splitLength, localLength);
			fullSplitReads.add(remote);
			assertEquals(S(r.getReadBases()), S(remote.getReadBases()));
		}
	}
	private SAMRecord sr(BreakpointSummary bs, int localLength, int splitLength) {
		SAMRecord r = Read(bs.referenceIndex, getLeftmostBasePos(bs, localLength, 0), bs.direction == FWD ? String.format("%dM%dS", localLength, splitLength) : String.format("%dS%dM", splitLength, localLength));
		String localBases = S(getRef(bs.referenceIndex, localLength, bs.start, bs.direction));
		String remoteBases = S(getRef(bs.remoteBreakend().referenceIndex, splitLength, bs.remoteBreakend().start, bs.remoteBreakend().direction));
		if (bs.direction == FWD) r.setReadBases(B(localBases + remoteBases)); 
		else r.setReadBases(B(remoteBases + localBases));
		return r;
	}
	/**
	 * Adds discordant read pairs to the evidence set
	 * @param reads output list
	 * @param dpExpectedReads 
	 * @param readLength
	 * @param offset offset from bp, zero indicates that the read aligns right up to the breakpoint
	 * @param fragSize
	 * @param fragSpan 
	 */
	private void addDP(List<SAMRecord> reads, List<SAMRecord> dpExpectedReads, BreakpointSummary bs, int readLength, int offset, int fragSize, int fragSpan) {
		assert(bs.referenceIndex == bs.referenceIndex2);
		int fragBasesAtLocal = readLength + offset;
		int fragBasesAtRemote = fragSize - fragBasesAtLocal;
		int unsequencedBasesAtRemote = fragBasesAtRemote - readLength;
		assert(fragBasesAtLocal >= readLength);
		assert(fragBasesAtRemote >= readLength);
		assert(bs.direction == FWD);
		SAMRecord[] dp = DP(bs.referenceIndex, getLeftmostBasePos(bs, readLength, offset), String.format("%dM", readLength), bs.direction == FWD,
				bs.remoteBreakend().referenceIndex, getLeftmostBasePos(bs.remoteBreakend(), fragBasesAtRemote, unsequencedBasesAtRemote), String.format("%dM", readLength), bs.remoteBreakend().direction == FWD);
		dp[0].setReadBases(B(S(RANDOM).substring(dp[0].getAlignmentStart() - 1,  dp[0].getAlignmentEnd())));
		dp[1].setReadBases(B(S(RANDOM).substring(dp[1].getAlignmentStart() - 1,  dp[1].getAlignmentEnd())));
		dp[0].setProperPairFlag(false);
		dp[1].setProperPairFlag(false);
		
		reads.add(dp[0]);
		reads.add(dp[1]);
		addExpectedPositions(dpExpectedReads, fragSize, fragSpan, dp[0], dp[1]);
		addExpectedPositions(dpExpectedReads, fragSize, fragSpan, dp[1], dp[0]);
	}
	public void addExpectedPositions(List<SAMRecord> dpExpectedReads, int fragSize, int fragSpan, SAMRecord anchor, SAMRecord read) {
		try {
			for (int i = -fragSpan; i <= fragSpan; i++) {
				int currentFragmentSize = fragSize - i;
				read = (SAMRecord)read.clone();
				boolean expectNegative = !anchor.getReadNegativeStrandFlag();
				if (read.getReadNegativeStrandFlag() != expectNegative) {
					read.setReadNegativeStrandFlag(expectNegative);
					read.setReadBases((byte[])read.getReadBases().clone());
					SequenceUtil.reverseComplement(read.getReadBases());
					read.setBaseQualities((byte[])read.getBaseQualities().clone());
					SequenceUtil.reverseComplement(read.getBaseQualities());
				}
				read.setReadUnmappedFlag(false);
				read.setProperPairFlag(true);
				read.setCigarString(String.format("%dM", read.getReadLength()));
				read.setReferenceIndex(anchor.getReferenceIndex());
				if (anchor.getReadNegativeStrandFlag()) {
					// expect read to be before
					read.setAlignmentStart(anchor.getAlignmentEnd() + SAMRecordUtil.getEndSoftClipLength(anchor) - currentFragmentSize + 1); 
				} else {
					read.setAlignmentStart(anchor.getAlignmentStart() - SAMRecordUtil.getStartSoftClipLength(anchor) + currentFragmentSize - read.getReadLength());
				}
				dpExpectedReads.add(read);
			}
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	/**
	 * Adds discordant read pairs to the evidence set
	 * @param reads output list
	 * @param readLength
	 * @param offset offset from bp, zero indicates that the read aligns right up to the breakpoint
	 * @param fragSize
	 */
	private void addRP(List<SAMRecord> reads, int readLength, int fragSize, int referenceIndex, int startPos) {
		SAMRecord[] rp = RP(referenceIndex, startPos, startPos + fragSize - readLength, readLength);
		rp[0].setReadBases(B(S(RANDOM).substring(rp[0].getAlignmentStart() - 1,  rp[0].getAlignmentEnd())));
		rp[1].setReadBases(B(S(RANDOM).substring(rp[1].getAlignmentStart() - 1,  rp[1].getAlignmentEnd())));
		reads.add(rp[0]);
		reads.add(rp[1]);
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
		SAMRecord[] dp = OEA(bs.referenceIndex, getLeftmostBasePos(bs, readLength, offset), String.format("%dM", readLength), bs.direction == FWD);
		dp[0].setReadBases(B(S(RANDOM).substring(dp[0].getAlignmentStart() - 1,  dp[0].getAlignmentEnd())));
		// TODO: reverseComp this?
		dp[1].setReadBases(B(S(getRef(bs.referenceIndex, localSplitBases, bs.start, bs.direction))
				+ S(getRef(bs.remoteBreakend().referenceIndex, fragBasesAtRemote, bs.remoteBreakend().start, bs.remoteBreakend().direction))));
		reads.add(dp[0]);
		reads.add(dp[1]);
	}
	//@Test
	public void generateMissassemblyGraph() throws IOException {
		// [1,10] GTAAAAAC
		//fail();
	}
	@Test
	public void small4merExample() throws IOException {
		BreakpointSummary bs = new BreakpointSummary(2, FWD, 100, 100, 2, BWD, 125, 125);
		int k = 4;
		int readLength = 10;
		int fragSize = 25;
		int fragSpan = 1;
		ProcessingContext pc = getCommandlineContext(false);
		
		pc.getConfig().getAssembly().k = k;
		pc.getConfig().getAssembly().writeFiltered = true;
		pc.getConfig().getAssembly().errorCorrection.maxBaseMismatchForCollapse = 0;
		pc.getConfig().getVariantCalling().writeFiltered = true;
		pc.getConfig().getVisualisation().directory = new File(super.testFolder.getRoot().getAbsoluteFile(), "visualisation");
		pc.getConfig().getVisualisation().assembly = true;
		pc.getConfig().getVisualisation().fullSizeAssembly = true;
		pc.getConfig().getVisualisation().assemblyProgress = true;
		pc.getConfig().getVisualisation().directory.mkdir();
				
		// Default 
		List<SAMRecord> variantReads = new ArrayList<SAMRecord>();
		List<SAMRecord> splitReadsRealignment = new ArrayList<SAMRecord>();
		List<SAMRecord> fullSplitReads = new ArrayList<SAMRecord>();
		List<SAMRecord> dpExpectedReads = new ArrayList<SAMRecord>();
		
		addDP(variantReads, dpExpectedReads, bs, readLength, 0, fragSize, fragSpan);
		addDP(variantReads, dpExpectedReads, bs, readLength, 4, fragSize - 1, fragSpan);
		
		addSR(variantReads, splitReadsRealignment, fullSplitReads, bs, readLength, 4, 0);
		addSR(variantReads, splitReadsRealignment, fullSplitReads, bs, readLength, 5, 0);
		addSR(variantReads, splitReadsRealignment, fullSplitReads, bs, readLength, 6, 0);
		// TODO: change soft clip position of read with sequencing error
		variantReads.get(variantReads.size()-1).getReadBases()[4] = 't';
		splitReadsRealignment.get(splitReadsRealignment.size()-1).getReadBases()[4] = 't';
		
		List<SAMRecord> baseReads = new ArrayList<SAMRecord>();
		Random random = new Random(0);
		for (int i = 1; i < 200; i++) {
			if (random.nextDouble() > 0.75) {
				addRP(baseReads, readLength, fragSize - fragSpan + random.nextInt(2 * fragSpan + 1), 2, i);
			}
		}
		
		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		reads.addAll(baseReads);
		reads.addAll(variantReads);
		createInput(reads);
		Collections.sort(splitReadsRealignment, new Ordering<SAMRecord>() {
			@Override
			public int compare(SAMRecord left, SAMRecord right) {
				return ComparisonChain.start()
				.compare(BreakpointFastqEncoding.getEncodedReferenceIndex(left.getReadName()), BreakpointFastqEncoding.getEncodedReferenceIndex(right.getReadName()))
				.compare(BreakpointFastqEncoding.getEncodedStartPosition(left.getReadName()), BreakpointFastqEncoding.getEncodedStartPosition(right.getReadName()))
				.result();
			}
		});
		// assemble evidence
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0, fragSize - 1, fragSize + 1);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(pc.getFileSystemContext().getRealignmentBam(input, 0), SortOrder.unsorted, splitReadsRealignment.toArray(new SAMRecord[0]));
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), output);
		aes.ensureAssembled();
		
		createBAM(new File(input.getParentFile(), "variant.bam"), SortOrder.coordinate, variantReads);
		createBAM(new File(input.getParentFile(), "dpexpected.bam"), SortOrder.coordinate, dpExpectedReads);
		createBAM(new File(input.getParentFile(), "srrealign.bam"), SortOrder.coordinate, splitReadsRealignment);
		createBAM(new File(input.getParentFile(), "srremote.bam"), SortOrder.coordinate, fullSplitReads);
		FileUtils.copyDirectory(input.getParentFile(), new File("C:/temp/tiny4"));
		System.err.append("tiny4merExample written to " + output.toString());
	}
}
