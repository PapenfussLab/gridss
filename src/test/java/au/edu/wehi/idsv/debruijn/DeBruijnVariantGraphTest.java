package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.EvidenceSubset;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SmallIndelSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnReadGraph;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class DeBruijnVariantGraphTest extends TestHelper {
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
		context.getAssemblyParameters().minReads = 0;
		context.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse = 0;
		context.getAssemblyParameters().errorCorrection.collapseBubblesOnly = true;
		context.getAssemblyParameters().writeFiltered = true;
		aes = AES(context);
		ass = new DeBruijnSubgraphAssembler(aes);
		result = Lists.newArrayList(); 
	}
	public void setup(int k, int maxPathTraversalNodes) {
		context = getContext();
		context.getAssemblyParameters().k = k;
		context.getAssemblyParameters().minReads = 0;
		context.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse = 0;
		context.getAssemblyParameters().errorCorrection.collapseBubblesOnly = true;
		context.getAssemblyParameters().writeFiltered = true;
		context.getAssemblyParameters().subgraph.traveralMaximumPathNodes = 0;
		aes = AES(context);
		ass = new DeBruijnSubgraphAssembler(aes);
		result = Lists.newArrayList(); 
	}
	public void setupWithGraphReduction(int k) {
		context = getContext();
		context.getAssemblyParameters().k = k;
		context.getAssemblyParameters().minReads = 0;
		context.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse = 1;
		context.getAssemblyParameters().errorCorrection.collapseBubblesOnly = false;
		context.getAssemblyParameters().writeFiltered = true;
		aes = AES(context);
		ass = new DeBruijnSubgraphAssembler(aes);
		result = Lists.newArrayList(); 
	}
	private static SAMRecord R(String read) {
		return R(null, read, null, false, true);
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
		SAMRecord local = Read(0, 1, "1M");
		remote.setReferenceIndex(1);
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
	public void should_assemble_single_dp_read() {
		addRead(R("AAAACGTC"), false);
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
	public void should_assemble_sc_with_dp() {
		addRead(R("CTAAA"), true);
		addRead(R(null, "AAAGT", null, false, true), false);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		
		assertEquals("CTAAAGT", S(result.get(0).getAssemblySequence()));
		assertEquals("AGT", S(result.get(0).getBreakendSequence()));
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
		
		SAMRecord[] oea = withSequence(unanchorSeq, RP(0, 1, 1000, unanchorSeq.length()));
		oea[0].setReadNegativeStrandFlag(mappedNegativeStrand);
		oea[1].setMateNegativeStrandFlag(mappedNegativeStrand);
		oea[1].setReadNegativeStrandFlag(mateNegativeStrand);
		oea[0].setMateNegativeStrandFlag(mateNegativeStrand);
		oea[0].setProperPairFlag(false);
		oea[1].setProperPairFlag(false);
		
		// Assembly should not depend on whether the read is mapped or not 
		setup(3);
		oea[0].setMateUnmappedFlag(true);
		oea[1].setReadUnmappedFlag(true);
		ass.addEvidence(SCE(direction, anchor));
		ass.addEvidence(NRRP(oea));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals(expectedSeq, S(result.get(0).getAssemblySequence()));
		assertEquals(breakendSequence, S(result.get(0).getBreakendSequence()));
		
		setup(3);
		oea[0].setMateUnmappedFlag(false);
		oea[1].setReadUnmappedFlag(false);
		ass.addEvidence(SCE(direction, anchor));
		ass.addEvidence(NRRP(oea));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals(expectedSeq, S(result.get(0).getAssemblySequence()));
		assertEquals(direction, result.get(0).getBreakendSummary().direction);
		assertEquals(breakendSequence, S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_assemble_end_soft_clipped_read() {
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		addRead(sc, true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("AAAACGTC", S(result.get(0).getAssemblySequence()));
		assertEquals("CGTC", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_assemble_start_soft_clipped_read() {
		SAMRecord sc = R("CGTCAAAA");
		sc.setCigarString("4S4M");
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, sc))));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("CGTCAAAA", S(result.get(0).getAssemblySequence()));
		assertEquals("CGTC", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void kmers_should_be_weighted_by_min_base_qual() {
		addRead(R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		// pad out read qualities
		assertArrayEquals(new byte[] { 4,5,5,5,6 }, ((SAMRecordAssemblyEvidence)result.get(0)).getSAMRecord().getBaseQualities());
	}
	@Test
	public void getBaseQuals_should_take_average_of_kmers_containing_base() {
		setup(3);
		SAMRecord sc = R(null, "CGTCAATTG", new byte[] {1,2,3,4,5,6,7,8,9}, false, true);
		sc.setCigarString("4M5S");
		addRead(sc, true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		// 123456789 k=3
		//  1234567
		assertArrayEquals(
				new byte[] {(1)/1,
						(1+2)/2,
						(1+2+3)/3,
						(2+3+4)/3,
						(3+4+5)/3,
						(4+5+6)/3,
						(5+6+7)/3,
						(6+7)/2,
						(7)/1,
				}, ((SAMRecordAssemblyEvidence)result.get(0)).getSAMRecord().getBaseQualities());
	}
	@Test
	public void should_greedy_traverse_highest_weight_path() {
		setup(3, 0);
		SAMRecord sc = R(null, "AAATCT", new byte[] { 1,2,3,4,5,6}, false, true);
		SAMRecord sc2 = R(null, "AAATCG", new byte[] { 2,3,4,5,6,7}, false, true);
		sc.setCigarString("3M3S");
		sc2.setCigarString("3M3S");
		addRead(sc, true);
		addRead(sc2, true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("AAATCG", S(result.get(0).getAssemblySequence()));
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
	@Ignore("Treated as unanchored")
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
	public void read_count_should_be_number_of_breakend_reads() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource normal = SES(false);
		MockSAMEvidenceSource tumour = SES(true);
		pc.getAssemblyParameters().k = 4;
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(aes);
		List<DirectedEvidence> e = Lists.newArrayList();
		e.add(SCE(FWD, normal, withSequence("TTCTT", Read(0, 9, "4M1S")))); // only on ref path
		e.add(SCE(FWD, normal, withSequence("TCTTA", Read(0, 10, "4M1S"))));
		e.add(SCE(FWD, normal, withSequence("TCTTAA", Read(0, 10, "4M2S"))));
		e.add(SCE(FWD, normal, withSequence("TCTTAAA", Read(0, 10, "4M3S"))));
		e.add(SCE(FWD, normal, withSequence("TCTTAAAC", Read(0, 10, "4M4S"))));
		e.add(SCE(FWD, tumour, withSequence("TCTTAAACG", Read(0, 10, "4M5S"))));
		e.add(SCE(FWD, tumour, withSequence("TCTTAAACGG", Read(0, 10, "4M6S"))));
		e.add(NRRP(withSequence(SequenceUtil.reverseComplement("AAACG"), OEA(0, 5, "5M", true))));	
		e.add(NRRP(tumour, withSequence(SequenceUtil.reverseComplement("CAAACG"), OEA(0, 6, "6M", true))));
		for (DirectedEvidence ev : e) assertFalse(ass.addEvidence(ev).iterator().hasNext());
		AssemblyEvidence result = ass.endOfEvidence().iterator().next().hydrateEvidenceSet(e).annotateAssembly();
		
		assertEquals("test assumes this contruction - did I get the test case wrong?", "TTCTTAAACGG", S(result.getAssemblySequence()));
		assertEquals(1, result.getAssemblySupportCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(1, result.getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR));
		assertEquals(2, result.getAssemblySupportCountReadPair(EvidenceSubset.ALL));
		assertEquals(4, result.getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(2, result.getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR));
		assertEquals(6, result.getAssemblySupportCountSoftClip(EvidenceSubset.ALL));
	}
	@Test
	public void should_assemble_spanning_breakend_with_insertion() {
		setup(25);
		SAMRecord f = Read(0, 1, "50M100S");
		SAMRecord b = Read(0, 100, "100S50M");
		f.setReadBases(B(S(RANDOM).substring(0, 50) + S(RANDOM).substring(500, 600)));
		b.setReadBases(B(S(RANDOM).substring(550, 650) + S(RANDOM).substring(100, 150)));
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, f))));
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, b))));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertTrue(result.get(0) instanceof SmallIndelSAMRecordAssemblyEvidence);
		SmallIndelSAMRecordAssemblyEvidence e = (SmallIndelSAMRecordAssemblyEvidence) result.get(0);
		assertEquals(S(RANDOM).substring(500, 650), e.getUntemplatedSequence());
		assertEquals(new BreakpointSummary(0, FWD, 50, 50, 0, BWD, 100, 100), e.getBreakendSummary());
	}
	@Test
	public void should_assemble_spanning_breakend() {
		setup(25);
		SAMRecord f = Read(0, 1, "50M50S");
		SAMRecord b = Read(0, 100, "50S50M");
		f.setReadBases(B(S(RANDOM).substring(0, 50) + S(RANDOM).substring(100, 150)));
		b.setReadBases(B(S(RANDOM).substring(0, 50) + S(RANDOM).substring(100, 150)));
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, f))));
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, b))));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertTrue(result.get(0) instanceof SmallIndelSAMRecordAssemblyEvidence);
		SmallIndelSAMRecordAssemblyEvidence e = (SmallIndelSAMRecordAssemblyEvidence) result.get(0);
		assertEquals("", e.getUntemplatedSequence());
		assertEquals(new BreakpointSummary(0, FWD, 50, 50, 0, BWD, 100, 100), e.getBreakendSummary());
		assertEquals("50M49D50M", e.getBackingRecord().getCigarString());
	}
	/**
	 * If the non-reference kmer path less than k-1 in length
	 * and is anchored on both sides, then the non-reference
	 * kmers are likely to be FPs as all bases have reference
	 * support 
	 */
	@Test
	public void should_assemble_even_if_length_less_than_k_misassemblies() {
		setup(4);
		// AAAACGTA
		// MMMMMSSS
		// SSSSMMMM
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAACGTA", Read(0, 10, "5M3S"))))));
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("AAAACGTA", Read(0, 10, "3S5M"))))));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
	}
	@Test
	public void should_assemble_non_reference_primary_kmer_as_reference_if_alt_ref_kmer_exists() {
		setupWithGraphReduction(4);
		addRead(withSequence("GCTAGG", Read(0, 1, "4M2S"))[0], true);
		addRead(withSequence("GCTAGG", Read(0, 1, "4M2S"))[0], true);
		addRead(withSequence("GCTATG", Read(0, 1, "5M1S"))[0], true);
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertEquals("GCTAGG", S(result.get(0).getAssemblySequence()));
		assertEquals("G", S(result.get(0).getBreakendSequence()));
	}
	@Test
	public void should_merge_kmers() {
		//          1         2         3         4
		// 1234567890123456789012345678901234567890123456789
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAA
		// ***************         ***************
		//        12345678         12345678
		// start anchor position = 8 -> 8 + k-1 = 15
		// end anchor position = 25
		aes = AES();
		aes.getContext().getConfig().getAssembly().k = 8;
		DeBruijnReadGraph g = new DeBruijnReadGraph(aes, 2, null);
		String seq = S(RANDOM).substring(0, 0+15) + S(RANDOM).substring(24, 24+15);
		SoftClipEvidence f = SCE(FWD, withSequence(seq, Read(2, 1, "15M15S")));
		SoftClipEvidence b = SCE(BWD, withSequence(seq, Read(2, 25, "15S15M")));
		g.addEvidence(f);
		g.addEvidence(b);
		
		assertEquals(30 - 8 + 1, g.size());
		//DeBruijnSubgraphNode node = g.getKmer(K(S(RANDOM).substring(0, 8)));		
	}
	@Test
	public void should_create_spanning_assembly() {
		// 1234567890123
		// GCTAGGTATCTCG
		// MMMM--
		// ---------MMMM
		setup(4);
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("GCTAGG", Read(0, 1, "4M2S"))[0]))));
		result.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("GCTAGGTATCTCG", Read(0, 10, "9S4M"))[0]))));
		result.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, result.size());
		assertTrue(result.get(0) instanceof SmallIndelSAMRecordAssemblyEvidence);
		SmallIndelSAMRecordAssemblyEvidence e = (SmallIndelSAMRecordAssemblyEvidence)result.get(0);
		assertEquals("GCTAGGTATCTCG", S(e.getAssemblySequence()));
		assertEquals("4M5I5D4M", e.getBackingRecord().getCigarString());
		assertEquals(1, e.getSAMRecord().getAlignmentStart());
		assertEquals(10, e.getRemoteSAMRecord().getAlignmentStart());
		assertEquals("5S4M", e.getRemoteSAMRecord().getCigarString());
	}
}