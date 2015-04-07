package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.List;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.AssemblyMethod;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.EvidenceSubset;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.anchored.DeBruijnAnchoredGraph;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;

import com.google.common.collect.Lists;

public class DeBruijnReadGraphTest extends TestHelper {
	private List<AssemblyEvidence> result;
	private DeBruijnSubgraphAssembler ass;
	private AssemblyEvidenceSource aes;
	private ProcessingContext context;
	@Before
	public void setup(){
		setup(4);
	}
	public void setup(int k) {
		context = getContext();
		context.getAssemblyParameters().k = k;
		context.getAssemblyParameters().method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		context.getAssemblyParameters().minReads = 0;
		context.getAssemblyParameters().maxBaseMismatchForCollapse = 0;
		context.getAssemblyParameters().collapseBubblesOnly = true;
		context.getAssemblyParameters().writeFilteredAssemblies = true;
		aes = AES(context);
		ass = new DeBruijnSubgraphAssembler(context, aes);
		result = Lists.newArrayList(); 
	}
	
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
		record.setMappingQuality(40);
		return record;
	}
	private SAMRecord inferLocal(BreakendDirection dir, SAMRecord remote) {
		SAMRecord local = Read(1, 1, "1M");
		local.setReadName(remote.getReadName());
		local.setReadPairedFlag(true);
		local.setReadNegativeStrandFlag((dir == BWD) == remote.getReadNegativeStrandFlag());
		local.setReadPairedFlag(false);
		local.setMappingQuality(40);
		local.setFirstOfPairFlag(true);
		remote.setFirstOfPairFlag(false);
		local.setSecondOfPairFlag(false);
		remote.setSecondOfPairFlag(true);
		local.setReadPairedFlag(true);
		remote.setReadPairedFlag(true);
		return local;
	}
	private void addRead(SAMRecord r, boolean sc) {
		if (sc) {
			result.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, r))));
		} else {
			result.addAll(Lists.newArrayList(ass.addEvidence(NRRP(SES(1000, 1000), inferLocal(FWD, r), r))));
		}
	}
	@Test
	public void should_assemble_single_read() {
		addRead(R("AAAACGTC"), true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		
		assertEquals("AAAACGTC", S(result.get(0).getAssemblySequence()));
	}
	@Test
	public void should_assemble_positive_strand_consensus() {
		addRead(R(null, "AAAACGTC", null, true, true), true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals(1, result.size());
		assertEquals("AAAACGTC", S(result.get(0).getAssemblySequence()));
	}
	@Test
	public void should_assemble_short_breakend() {
		addRead(withSequence("CTAAA", Read(0, 1, "4M1S"))[0], true);
		addRead(withSequence("CTAAA", Read(0, 2, "4M1S"))[0], true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		
		assertEquals("CTAAA", S(result.get(0).getAssemblySequence()));
		assertEquals("A", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_assemble_unanchored_reads() {
		addRead(withSequence("CTAAA", Read(0, 1, "4M1S"))[0], true);
		addRead(R(null, "AAAGT", null, false, true), false);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		
		assertEquals("CTAAAGT", S(result.get(0).getAssemblySequence()));
		assertEquals("AGT", S(result.get(0).getBreakendSequence()));
	}
	@Test  
	public void unanchored_reads_should_require_mapped_mate() {
		addRead(R("CTAAA"), true);
		SAMRecord unanchored = R(null, "AAAGT", null, false, true);
		unanchored.setMateUnmappedFlag(true);
		addRead(unanchored, false);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		
		assertEquals("CTAAAGT", S(result.get(0).getAssemblySequence()));
		assertEquals("AAAGT", S(result.get(0).getBreakendSequence()));
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
		setup(3);
		addRead(anchor, true);
		SAMRecord unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(true);
		addRead(unanchored, false);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals(expectedSeq, S(result.get(0).getAssemblySequence()));
		assertEquals(breakendSequence, S(result.get(0).getBreakendSequence()));
		
		setup(3);
		addRead(anchor, true);
		unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(false);
		addRead(unanchored, false);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals(expectedSeq, S(result.get(0).getAssemblySequence()));
		assertEquals(direction, result.get(0).getBreakendSummary().direction);
		assertEquals(breakendSequence, S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_assemble_soft_clipped_read() {
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");e
		addRead(sc, true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("AAAACGTC", S(result.get(0).getAssemblySequence()));
		assertEquals("CGTC", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_greedy_traverse_highest_weight_path() {
		setup(2);
		SAMRecord sc = R(null, "ACGTACTGAG", new byte[] { 1,2,3,4,5,6,7,8,9,10}, false, true);
		sc.setCigarString("4M6S");
		addRead(sc, true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("TACTGAGT", S(result.get(0).getAssemblySequence()));
	}
	@Test
	public void should_assemble_backward_breakpoint() {
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4S4M");
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, sc))));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("AAAACGTC", S(result.get(0).getAssemblySequence()));
		assertEquals("AAAA", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_use_offset_kmer_if_softclip_longer_than_k() {
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("2M6S");
		addRead(sc, true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("AAAACGTC", S(result.get(0).getAssemblySequence()));
		assertEquals("AACGTC", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void assembly_base_quality_should_be_sum_of_min_read_qualities() {
		addRead(R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		// pad out read qualities
		assertArrayEquals(new byte[] { /*4,6,8,8,*/8 }, result.get(0).getBreakendQuality());
	}
	@Test // can't see the anchor qualities in the current API
	public void assembly_base_quality_should_pad_to_match_read_length() {
		addRead(R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		// pad out read qualities
		assertEquals(S(result.get(0).getAssemblySequence()), result.get(0).getBreakendQuality());
	}
	@Test
	public void read_base_count_should_be_number_of_breakend_read_bases() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource normal = SES(false);
		MockSAMEvidenceSource tumour = SES(true);
		pc.getAssemblyParameters().k = 4;
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(pc, aes);
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTA", Read(0, 10, "4M1S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTAA", Read(0, 10, "4M2S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTAAA", Read(0, 10, "4M3S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTTTAAAC", Read(0, 10, "4M4S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TTTTAAACG", Read(0, 10, "4M5S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TTTTAAACGG", Read(0, 10, "4M6S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(NRRP(withSequence("CGTTT", OEA(0, 5, "5M", true)))).iterator().hasNext());	
		assertFalse(ass.addEvidence(NRRP(tumour, withSequence("CGTTTG", OEA(0, 6, "5M", true)))).iterator().hasNext());
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
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(pc, aes);
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TTCTT", Read(0, 9, "4M1S")))).iterator().hasNext()); // only on ref path
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTA", Read(0, 10, "4M1S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTAA", Read(0, 10, "4M2S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTAAA", Read(0, 10, "4M3S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, normal, withSequence("TCTTAAAC", Read(0, 10, "4M4S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TCTTAAACG", Read(0, 10, "4M5S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(SCE(FWD, tumour, withSequence("TCTTAAACGG", Read(0, 10, "4M6S")))).iterator().hasNext());
		assertFalse(ass.addEvidence(NRRP(withSequence(SequenceUtil.reverseComplement("AAACG"), OEA(0, 5, "5M", true)))).iterator().hasNext());	
		assertFalse(ass.addEvidence(NRRP(tumour, withSequence(SequenceUtil.reverseComplement("CAAACG"), OEA(0, 6, "5M", true)))).iterator().hasNext());
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