package au.edu.wehi.idsv.debruijn.anchoured;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.EvidenceSubset;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.anchored.DeBruijnAnchoredGraph;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;

public class DeBruijnReadGraphTest extends TestHelper {
	private static SAMRecord R(String read) {
		return R(null, read, null, false, true);
	}
	private static SAMRecord R(String readName, String read) {
		return R(readName, read, null, false, true);
	}
	private static SAMRecord R(String readName, String read, byte[] qual, boolean mappedNegativeStrand, boolean mateNegativeStrand) {
		SAMRecord record = new SAMRecord(getHeader());
		if (qual == null) {
			qual = new byte[read.length()];
			for (int i = 0; i < qual.length; i++) qual[i] = 1;
		}
		record.setReadBases(B(read));
		record.setBaseQualities(qual);
		record.setReadPairedFlag(true);
		record.setReadNegativeStrandFlag(mappedNegativeStrand);
		record.setReadUnmappedFlag(false);
		record.setMateUnmappedFlag(false);
		record.setMateNegativeStrandFlag(mateNegativeStrand);
		if (readName == null) {
			readName = String.format("%s-%s-%s%s", read, qual, mappedNegativeStrand, mateNegativeStrand);
		}
		record.setReadName(readName);
		record.setReferenceIndex(0);
		record.setAlignmentStart(1);
		record.setCigarString(String.format("%dM1S", read.length() - 1));
		return record;
	}
	private SAMRecord inferLocal(DeBruijnAnchoredGraph ass, SAMRecord remote) {
		SAMRecord local = Read(0, 1, "1M");
		local.setReadName(remote.getReadName());
		local.setReadPairedFlag(true);
		local.setReadNegativeStrandFlag(ass.getDirection() == BWD);
		return local;
	}
	private void addRead(DeBruijnAnchoredGraph ass, SAMRecord r, boolean sc) {
		if (sc) {
			ass.addEvidence(SCE(ass.getDirection(), r));
		} else {
			ass.addEvidence(NRRP(inferLocal(ass, r), r));
		}
	}
	private void removeRead(DeBruijnAnchoredGraph ass, SAMRecord r, boolean sc) {
		if (sc) {
			ass.removeEvidence(SCE(ass.getDirection(), r));
		} else {
			ass.removeEvidence(NRRP(inferLocal(ass, r), r));
		}
	}
	@Test
	public void should_assemble_single_read() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(), 4, BreakendDirection.Forward);
		addRead(ass, R("AAAACGTC"), true);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("AAAACGTC", S(result.getAssemblySequence()));
	}
	@Test
	public void should_assemble_positive_strand_consensus() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(), 4, BreakendDirection.Forward);
		addRead(ass, R(null, "AAAACGTC", null, true, true), true);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("AAAACGTC", S(result.getAssemblySequence()));
	}
	@Test
	public void should_assemble_unanchored_reads() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(),3, BreakendDirection.Forward);
		addRead(ass, withSequence("CTAAA", Read(0, 1, "4M1S"))[0], true);
		addRead(ass, R(null, "AAAGT", null, false, true), false);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("CTAAAGT", S(result.getAssemblySequence()));
		assertEquals("AGT", S(result.getBreakendSequence()));
	}
	@Test(expected = RuntimeException.class)  
	public void unanchored_reads_should_require_mapped_mate() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(),3, BreakendDirection.Forward);
		addRead(ass, R("CTAAA"), true);
		SAMRecord unanchored = R(null, "AAAGT", null, false, true);
		unanchored.setMateUnmappedFlag(true);
		addRead(ass, unanchored, false);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("CTAAAGT", S(result.getAssemblySequence()));
		assertEquals("AAAGT", S(result.getBreakendSequence()));
	}
	@Test
	public void should_assemble_unanchored_reads_in_FR_orientation() {
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "AGT", "CTAAA", "AAAGT", true, false);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "AGT", "CTAAA", "AAAGT", false, true);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "AGT", "CTAAA", "ACTTT", true, true);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "AGT", "CTAAA", "ACTTT", false, false);
		
		assertExpected(BreakendDirection.Backward, "GTAAACT", "GTA", "AAACT", "GTAAA", true, false);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "GTA", "AAACT", "GTAAA", false, true);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "GTA", "AAACT", "TTTAC", true, true);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "GTA", "AAACT", "TTTAC", false, false);
	}
	private void assertExpected(BreakendDirection direction, String expectedSeq, String breakendSequence , String anchorSeq, String unanchorSeq, boolean mappedNegativeStrand, boolean mateNegativeStrand) {
		SAMRecord anchor = R(anchorSeq);
		if (direction == BWD) anchor.setCigarString("1S4M");
		
		// Assembly should not depend on whether the read is mapped or not 
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(),3, direction);
		addRead(ass, anchor, true);
		SAMRecord unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(true);
		addRead(ass, unanchored, false);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals(expectedSeq, S(result.getAssemblySequence()));
		assertEquals(breakendSequence, S(result.getBreakendSequence()));
		
		ass = new DeBruijnAnchoredGraph(getContext(), AES(),3, direction);
		addRead(ass, anchor, true);
		unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(false);
		addRead(ass, unanchored, false);
		result = ass.assemble(0, 1);
		assertEquals(expectedSeq, S(result.getAssemblySequence()));
		assertEquals(direction, result.getBreakendSummary().direction);
		assertEquals(breakendSequence, S(result.getBreakendSequence()));
	}
	@Test
	public void should_assemble_soft_clipped_read() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(), 4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		addRead(ass, sc, true);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("AAAACGTC", S(result.getAssemblySequence()));
		assertEquals("CGTC", S(result.getBreakendSequence()));
	}
	@Test
	public void should_greedy_traverse_highest_weight_path() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(),2, BreakendDirection.Forward);
		SAMRecord sc = R(null, "ACGTACTGAG", new byte[] { 1,2,3,4,5,6,7,8,9,10}, false, true);
		sc.setCigarString("4M6S");
		addRead(ass, sc, true);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("TACTGAGT", S(result.getAssemblySequence()));
	}
	@Test
	public void should_assemble_backward_breakpoint() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(), 4, BreakendDirection.Backward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4S4M");
		addRead(ass, sc, true);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("AAAACGTC", S(result.getAssemblySequence()));
		assertEquals("AAAA", S(result.getBreakendSequence()));
	}
	@Test
	public void remove_should_exclude_from_graph() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(), 4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		addRead(ass, sc, true);
		addRead(ass, R("r1", "AAAATTTT"), false);
		addRead(ass, R("r2", "AAAATTTT"), false);
		assertEquals("AAAATTTT", S(ass.assemble(0, 1).getAssemblySequence()));
		removeRead(ass, R("r1", "AAAATTTT"), false);
		removeRead(ass, R("r2", "AAAATTTT"), false);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("AAAACGTC", S(result.getAssemblySequence()));
	}
	@Test
	public void should_use_offset_kmer_if_softclip_longer_than_k() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(), 4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("2M6S");
		addRead(ass, sc, true);
		AssemblyEvidence result = ass.assemble(0, 1);
		assertEquals("AAAACGTC", S(result.getAssemblySequence()));
		assertEquals("AACGTC", S(result.getBreakendSequence()));
	}
	@Test
	public void assembly_base_quality_should_be_sum_of_min_read_qualities() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(),3, BreakendDirection.Forward);
		addRead(ass, R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(ass, R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		AssemblyEvidence result = ass.assemble(0, 1);
		// pad out read qualities
		assertArrayEquals(new byte[] { /*4,6,8,8,*/8 }, result.getBreakendQuality());
	}
	// @Test // can't see the anchor qualities in the current API
	public void assembly_base_quality_should_pad_to_match_read_length() {
		DeBruijnAnchoredGraph ass = new DeBruijnAnchoredGraph(getContext(), AES(),3, BreakendDirection.Forward);
		addRead(ass, R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(ass, R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		AssemblyEvidence result = ass.assemble(0, 1);
		// pad out read qualities
		assertEquals(S(result.getAssemblySequence()), result.getBreakendQuality());
	}
	@Test
	public void read_base_count_should_be_number_of_breakend_read_bases() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource normal = SES(false);
		MockSAMEvidenceSource tumour = SES(true);
		pc.getAssemblyParameters().k = 4;
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(pc, AES());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTA", Read(0, 10, "4M1S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTAA", Read(0, 10, "4M2S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTAAA", Read(0, 10, "4M3S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTAAAC", Read(0, 10, "4M4S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TTTTAAACG", Read(0, 10, "4M5S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TTTTAAACGG", Read(0, 10, "4M6S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(NRRP(withSequence("CGTTT", OEA(0, 5, "5M", true)))).iterator().hasNext());	
		assertFalse(ass.addEvidence(NRRP(tumour, withSequence("CGTTTG", OEA(0, 5, "5M", true)))).iterator().hasNext());
		AssemblyEvidence result = ass.endOfEvidence().iterator().next();
		
		assertEquals("test assumes this contruction - did I get the test case wrong?", "TTTTAAACGG", S(result.getAssemblySequence()));
		assertEquals(5 + 6 + 5, result.getAssemblyBaseCount(EvidenceSubset.TUMOUR));
		assertEquals(1 + 2 + 3 + 4 + 5, result.getAssemblyBaseCount(EvidenceSubset.NORMAL));
	}
	@Test
	public void read_count_should_be_number_of_breakend_reads() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource normal = SES(false);
		MockSAMEvidenceSource tumour = SES(true);
		pc.getAssemblyParameters().k = 4;
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(pc, AES());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTCTT", Read(0, 9, "4M1S")))).iterator().hasNext()); // only on ref path
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTA", Read(0, 10, "4M1S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTAA", Read(0, 10, "4M2S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTAAA", Read(0, 10, "4M3S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTAAAC", Read(0, 10, "4M4S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TCTTAAACG", Read(0, 10, "4M5S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TCTTAAACGG", Read(0, 10, "4M6S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(NRRP(withSequence(SequenceUtil.reverseComplement("AAACG"), OEA(0, 5, "5M", true)))).iterator().hasNext());	
		assertFalse(ass.addEvidence(NRRP(tumour, withSequence(SequenceUtil.reverseComplement("CAAACG"), OEA(0, 5, "5M", true)))).iterator().hasNext());
		AssemblyEvidence result = ass.endOfEvidence().iterator().next();
		
		assertEquals("test assumes this contruction - did I get the test case wrong?", "TTCTTAAACGG", S(result.getAssemblySequence()));
		assertEquals(1, result.getAssemblySupportCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(1, result.getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR));
		assertEquals(2, result.getAssemblySupportCountReadPair(EvidenceSubset.ALL));
		assertEquals(4, result.getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(2, result.getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR));
		assertEquals(6, result.getAssemblySupportCountSoftClip(EvidenceSubset.ALL));
	}
}