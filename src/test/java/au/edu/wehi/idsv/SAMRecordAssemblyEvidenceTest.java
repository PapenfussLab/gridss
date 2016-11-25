package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.configuration.ConfigurationException;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;


public class SAMRecordAssemblyEvidenceTest extends TestHelper {
	@Test
	public void should_create_SAMRecord_for_assembly() {
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		assertNotNull(e);
		assertEquals("GTAC", S(e.getReadBases()));
		assertArrayEquals( new byte[] {1,2,3,4}, e.getBaseQualities());
	}
	@Test
	public void anchor_positions_should_match_genomic() {
		SAMRecord fwd = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null, 1, 10, 3, B("GTACCCA"), new byte[] { 1, 2, 3, 4, 4, 8, 8 });
		SAMRecord bwd = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null, 1, 10, 3, B("GTACCCA"), new byte[] { 1, 2, 3, 4, 4, 8, 8 });
		assertEquals(1, (int)fwd.getReferenceIndex());
		assertEquals(1, (int)bwd.getReferenceIndex());
		assertEquals(8, fwd.getAlignmentStart());
		assertEquals(10, bwd.getAlignmentStart());
		assertEquals("3M4S", fwd.getCigarString());
		assertEquals("4S3M", bwd.getCigarString());
		assertEquals("CCCA", S(asEvidence(fwd).getBreakendSequence()));
		assertEquals("GTAC", S(asEvidence(bwd).getBreakendSequence()));
		assertEquals("GTA", S(asEvidence(fwd).getAnchorSequence()));
		assertEquals("CCA", S(asEvidence(bwd).getAnchorSequence()));
	}
	@Test
	public void unanchor_positions_should_match_genomic() {
		SAMRecord fwd = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 7, 5, 10), null, B("AAA"), B("AAA"), new int[] {2, 0});
		SAMRecord bwd = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, BWD, 7, 5, 10), null, B("AAA"), B("AAA"), new int[] {2, 0});
		assertEquals(1, (int)fwd.getReferenceIndex());
		assertEquals(1, (int)bwd.getReferenceIndex());
		assertEquals(5, fwd.getAlignmentStart());
		assertEquals(5, bwd.getAlignmentStart());
		assertEquals(10, fwd.getAlignmentEnd());
		assertEquals(10, bwd.getAlignmentEnd());
		assertEquals("1X4N1X3S", fwd.getCigarString());
		assertEquals("3S1X4N1X", bwd.getCigarString());
		assertEquals("AAA", S(asEvidence(fwd).getBreakendSequence()));
		assertEquals("AAA", S(asEvidence(bwd).getBreakendSequence()));
		assertEquals("", S(asEvidence(fwd).getAnchorSequence()));
		assertEquals("", S(asEvidence(bwd).getAnchorSequence()));
	}
	@Test
	public void unanchor_sequences_should_match_assembly() {
		for (SAMRecord e : new SAMRecord[] {
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 1, 1, 1), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 1, 1, 2), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 1, 1, 3), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 1, 1, 4), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, BWD, 1, 1, 1), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, BWD, 1, 1, 2), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, BWD, 1, 1, 3), null, B("AAA"), B("AAA"), new int[] {2, 0}),
				AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, BWD, 1, 1, 4), null, B("AAA"), B("AAA"), new int[] {2, 0})
		}) {
			assertEquals("AAA", S(asEvidence(e).getBreakendSequence()));
			assertEquals("", S(asEvidence(e).getAnchorSequence()));
		}
	}
	private void assertPairing(SAMRecord assembly, SAMRecord realign) {
		assertNotNull(assembly);
		assertNotNull(realign);
		assertTrue(assembly.getReadPairedFlag());
		assertTrue(realign.getReadPairedFlag());
		assertTrue(assembly.getFirstOfPairFlag());
		assertFalse(realign.getFirstOfPairFlag());
		assertFalse(assembly.getSecondOfPairFlag());
		assertTrue(realign.getSecondOfPairFlag());
		assertEquals(assembly.getReadName(), realign.getReadName());
		assertMateFields(assembly, realign);
		assertMateFields(realign, assembly);
	}
	private void assertMateFields(SAMRecord r, SAMRecord mate) {
		assertEquals(r.getMateNegativeStrandFlag(), mate.getReadNegativeStrandFlag());
		assertEquals(r.getMateUnmappedFlag(), mate.getReadUnmappedFlag());
		assertEquals(r.getMateReferenceIndex(), mate.getReferenceIndex());
		assertEquals(r.getMateAlignmentStart(), mate.getAlignmentStart());
	}
	@Test
	public void should_use_MS_for_anchored_breakend() {
		assertEquals("1M3S", AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}).getCigarString());
		assertEquals("3S1M", AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null,
				1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}).getCigarString());
	}
	@Test
	public void should_use_XNXS_for_unanchored_interval_breakend() {
		assertEquals("1X1N1X4S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, FWD, 1, 1, 3), null, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0}).getCigarString());
		assertEquals("4S1X9N1X", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, BWD, 10, 10, 20), null, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0}).getCigarString());
	}
	@Test
	public void should_use_minimal_cigar_representation() {
		assertEquals("1X3S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 5, 5, 5), null, B("AAA"), B("AAA"), new int[] {2, 0}).getCigarString());
		assertEquals("2X3S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 5, 5, 6), null, B("AAA"), B("AAA"), new int[] {2, 0}).getCigarString());
		assertEquals("1X1N1X3S", AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(1, FWD, 5, 5, 7), null, B("AAA"), B("AAA"), new int[] {2, 0}).getCigarString());
	}
	public static void assertEvidenceEquals(SAMRecord ass1, SAMRecord ass2) {
		SingleReadEvidence e = asEvidence(ass1);
		SingleReadEvidence r = asEvidence(ass2);
		AssemblyAttributes ea = new AssemblyAttributes(e);
		AssemblyAttributes ra = new AssemblyAttributes(e);
		assertEquals(e.getAnchorSequence().length, r.getAnchorSequence().length);
		assertEquals(S(e.getAnchorSequence()) , S(r.getAnchorSequence()));
		assertEquals(ea.getAssemblyReadPairLengthMax(0) , ra.getAssemblyReadPairLengthMax(0));
		assertEquals(ea.getAssemblyReadPairLengthMax(1) , ra.getAssemblyReadPairLengthMax(1));
		assertEquals(ea.getAssemblyReadPairLengthMax() , ra.getAssemblyReadPairLengthMax());
		assertEquals(S(ass1.getBaseQualities()) , S(ass2.getBaseQualities()));
		assertEquals(ea.getAssemblySoftClipLengthMax(0) , ra.getAssemblySoftClipLengthMax(0));
		assertEquals(ea.getAssemblySoftClipLengthMax(1) , ra.getAssemblySoftClipLengthMax(1));
		assertEquals(ea.getAssemblySoftClipLengthMax() , ra.getAssemblySoftClipLengthMax());
		assertEquals(ea.getAssemblySoftClipLengthTotal(0) , ra.getAssemblySoftClipLengthTotal(0));
		assertEquals(ea.getAssemblySoftClipLengthTotal(1) , ra.getAssemblySoftClipLengthTotal(1));
		assertEquals(ea.getAssemblySoftClipLengthTotal() , ra.getAssemblySoftClipLengthTotal());
		assertEquals(ea.getAssemblySupportCountReadPair(0) , ra.getAssemblySupportCountReadPair(0));
		assertEquals(ea.getAssemblySupportCountReadPair(1) , ra.getAssemblySupportCountReadPair(1));
		assertEquals(ea.getAssemblySupportCountReadPair() , ra.getAssemblySupportCountReadPair());
		assertEquals(ea.getAssemblySupportCountSoftClip(0) , ra.getAssemblySupportCountSoftClip(0));
		assertEquals(ea.getAssemblySupportCountSoftClip(1) , ra.getAssemblySupportCountSoftClip(1));
		assertEquals(ea.getAssemblySupportCountSoftClip() , ra.getAssemblySupportCountSoftClip());
		assertArrayEquals(e.getBreakendQuality() , r.getBreakendQuality());
		assertEquals(S(e.getBreakendSequence()) , S(r.getBreakendSequence()));
		assertEquals(e.getBreakendSummary(), r.getBreakendSummary());
		if (e.getBreakendSummary() != null) {
			assertEquals(e.getBreakendSummary().getClass(), r.getBreakendSummary().getClass());
		}
		assertEquals(e.getEvidenceID() , r.getEvidenceID());
		assertEquals(e.getEvidenceSource() , r.getEvidenceSource());
		assertEquals(e.getLocalMapq() , r.getLocalMapq());
		if (e instanceof DirectedBreakpoint) {
			assertTrue(r instanceof DirectedBreakpoint);
			DirectedBreakpoint de = (DirectedBreakpoint)e;
			DirectedBreakpoint dr = (DirectedBreakpoint)r;
			assertEquals(de.getBreakendSummary(), dr.getBreakendSummary());
			assertEquals(de.getRemoteMapq(), dr.getRemoteMapq());
			assertEquals(de.getUntemplatedSequence(), dr.getUntemplatedSequence());
		}
	}
	@Test
	public void should_track_breakend_evidence() {
		DirectedEvidence e1 = SCE(BWD, Read(0, 1, "5S5M"));
		DirectedEvidence e2 = SCE(BWD, Read(0, 1, "6S5M"));
		DirectedEvidence e3 = NRRP(OEA(0, 1, "1M", false));
		DirectedEvidence b1 = SCE(BWD, Read(0, 1, "7S5M"));
		SAMRecord r = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.newArrayList(e1, e2, e3),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		AssemblyAttributes e = new AssemblyAttributes(r);
		assertTrue(e.isPartOfAssembly(e1));
		assertTrue(e.isPartOfAssembly(e2));
		assertTrue(e.isPartOfAssembly(e3));
		assertFalse(e.isPartOfAssembly(b1));
		
		assertEquals(2, e.getAssemblySupportCountSoftClip());
		assertEquals(1, e.getAssemblySupportCountReadPair());
	}
	@Test
	public void getEvidenceIDs_should_return_underlying_evidence() {
		DirectedEvidence e1 = SCE(BWD, Read(0, 1, "5S5M"));
		DirectedEvidence e2 = SCE(BWD, Read(0, 1, "6S5M"));
		DirectedEvidence e3 = NRRP(OEA(0, 1, "1M", false));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.newArrayList(e1, e2, e3),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		//Collection<String> ids = e.getEvidenceIDs();
		Collection<String> eid = new AssemblyAttributes(e).getEvidenceIDs();
		assertEquals(3, eid.size());
		assertTrue(eid.contains(e1.getEvidenceID()));
		assertTrue(eid.contains(e2.getEvidenceID()));
		assertTrue(eid.contains(e3.getEvidenceID()));
	}
	@Test
	public void getEvidenceIDs_should_return_empty_collection_for_no_evidence() {
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.newArrayList(),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4});
		//Collection<String> ids = e.getEvidenceIDs();
		assertEquals(0, new AssemblyAttributes(e).getEvidenceIDs().size());
	}
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	@Test
	public void realign_should_fix_778_chr1_170849702_misalignment() throws ConfigurationException {
		String assembly = "ATCCATCCCTATGACCCAAACATCTCCCACCAGGCCTCATGTTCAATATTAAAGATCACATTTCAACTTGAGATTTGGAGGGGACAAACATACAAATCATATCATTATCTCTCTCCCCACTTCTCTCTTTATCAATCCCTCCCTCTTTGTCAATCTTAGCCTTGGCCTTCAGATTTTACCACTTGATTTTTCACATTTTCTGTATTCTTAAT"
				+ "GATTATTATATTTTCATGTTCTTGCTAATCTATATCATGGTTAGAAATCAAAGCATGCCGAAATTTCTCTCTTACTTTTTTTGCTGTT";
		File ref = new File("src/test/resources/chr1_170849600_170849850.fa");
		ProcessingContext context = new ProcessingContext(getFSContext(), ref, null, new ArrayList<Header>(),
				new GridssConfiguration());
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(context, AES(context), BWD, null,
				0, 170849702-170849600+1, 97, B(assembly), B(40, assembly.length()));
		e = SAMRecordUtil.realign(context.getReference(), e, 50, true);
		// anchor location is 11bp off
		assertEquals("212S88M", e.getCigarString());
		assertEquals(170849713-170849600+1, asEvidence(e).getBreakendSummary().start);
	}
	@Test
	public void small_indel_should_be_called_if_realignment_spans_event() {
		String assembly = "AAAAAAAAAATTAAAAAAAAAA";
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				0, 1, 10, B(assembly), B(40, assembly.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("10M2I10M", e.getCigarString());
	}
	/**
	 * Needed for realignment that turn out to match reference
	 */
	@Test
	public void should_allow_reference_allele_assemblies() {
		String assembly = "AAAAAAAAAAA";
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				0, 1, assembly.length(), B(assembly), B(40, assembly.length()));
		
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), e).size());
	}
	@Test
	public void should_allow_realignment_to_reference_allele() {
		String assembly = "AAAAAAAAAAA";
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				0, 1, 1, B(assembly), B(40, assembly.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), e).size());
	}
	@Test
	public void getBreakendQual_should_exclude_assembled_evidence_that_does_not_support_breakend() {
		List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				NRRP(OEA(0, 1, "1M", true)),
				NRRP(OEA(0, 2, "1M", true)),
				NRRP(OEA(0, 1000, "1M", true)));
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().excludeNonSupportingEvidence = false;
		SAMRecord ass = AssemblyFactory.createUnanchoredBreakend(pc, AES(pc), new BreakendSummary(0, FWD, 1), support, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0});
		assertEquals(support.get(0).getBreakendQual() + support.get(1).getBreakendQual() + support.get(2).getBreakendQual(), asEvidence(ass).getBreakendQual(), DELTA);
		
		pc.getAssemblyParameters().excludeNonSupportingEvidence = true;
		ass = AssemblyFactory.createUnanchoredBreakend(pc, AES(pc), new BreakendSummary(0, FWD, 1), support, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0});
		assertEquals(support.get(0).getBreakendQual() + support.get(1).getBreakendQual(), asEvidence(ass).getBreakendQual(), DELTA);
	}
	@Test
	public void getAnchor_should_not_include_inexact_breakend_bases() {
		//List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(NRRP(OEA(0, 1, "1M", true)));
		assertEquals(0, asEvidence(AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new BreakendSummary(0, FWD, 1), null, B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0})).getAnchorSequence().length);
	}
	@Test
	public void breakpoint_should_use_indel_cigar() {
		//List<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList();
		// 1234567890   1234567890
		//         MMIIIDDDDMMMM
		//         NNAAA    TTTT
		SAMRecord r = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), null,
				0, 10, 2,
				0, 15, 4,
				B("NNAAATTTT"),
				B("ABCDEFGHI"));
		assertEquals(0, (int)r.getReferenceIndex());
		assertEquals(9, r.getAlignmentStart());
		assertEquals("2M3I4D4M", r.getCigarString());
		assertEquals("ABCDEFGHI", S(r.getBaseQualities()));
		assertEquals("NNAAATTTT", S(r.getReadBases()));
	}
	/*// Now covered by the more general SA split read realignment tests
	@Test
	public void getAllRealignments_should_return_all_breakpoints_fwd() {
		//                1         2         3         4
		//      012345678901234567890123456789012345678901234567890123456789012345
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA anchor
		//      ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA realign 1
		//      SSSSSSSSSSSSSSSSSSSSSSSSSMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
		//          CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT realign 2
		//          SSSSSSSSMMMMMMMMSSSSSSSSSSSSSSSSS
		//                                        G realign 3
		//                                        M
		//
		//
		//          1         2         3         4         5         6         7      
		// 123456789012345678901234567890123456789012345678901234567890123456789012345
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA
		// *****SSSSSSSSSSSSSSSSSSSSSSSSSMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
		//          SSSSSSSSMMMMMMMMSSSSSSSSSSSSSSSSS
		//                                        M
		SAMRecord be = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				0, 5, 5, B("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA"), B(40,75));
		RealignedSAMRecord e = (RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(getContext(), be, ImmutableList.of(
				withReadName("0#0#0#readname", withSequence(B("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA"),
						withQual(B(40,"ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA".length()), Read(1, 100, "25S5M40S"))))[0],
				withReadName("0#0#4#readname", withSequence(B("CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT"),
						withQual(B(40,"CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT".length()), Read(2, 200, "8S8M17S"))))[0],
				withReadName("0#0#34#readname", withQual(B("1"), withSequence("G", Read(0, 40, "1M"))))[0]
				));
		assertEquals(new BreakpointSummary(0, FWD, 5, 2, BWD, 200), e.getBreakendSummary());
		assertEquals("ATCGCAAGAGCG", e.getUntemplatedSequence());
		List<SAMRecordAssemblyEvidence> rl = e.getSubsequentRealignments();
		assertEquals(2, rl.size());
		
		RealignedSAMRecord r0 = (RealignedSAMRecordAssemblyEvidence)rl.get(0);
		assertEquals(new BreakpointSummary(2, FWD, 207, 1, BWD, 100), r0.getBreakendSummary());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTAT", S(r0.getAnchorSequence()));
		assertEquals("TCGAC", r0.getUntemplatedSequence());
		
		RealignedSAMRecord r1 = (RealignedSAMRecordAssemblyEvidence)rl.get(1);
		assertEquals(new BreakpointSummary(1, FWD, 104, 0, BWD, 40), r1.getBreakendSummary());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAA", S(r1.getAnchorSequence()));
		assertEquals("GTCA", r1.getUntemplatedSequence());
	}
	@Test
	public void getAllRealignments_should_return_all_breakpoints_fwd_fwd_bwd() {
		//          1         2         3         4         5         6         7      
		// 123456789012345678901234567890123456789012345678901234567890123456789012345
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA
		// *****SSSSSSSSSSSSSSSSSSSSSSSSSMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
		//          SSSSSSSSMMMMMMMMSSSSSSSSSSSSSSSSS
		//                                        M
		SAMRecord be = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				0, 5, 5, B("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA"), B(40,75));
		RealignedSAMRecord e = (RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(getContext(), be, ImmutableList.of(
				withReadName("0#0#0#readname", withSequence(B("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA"),
						withQual(B(40,"ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA".length()), Read(1, 100, "25S5M40S"))))[0],
				withReadName("0#0#4#readname", withSequence(B("CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT"),
						withQual(B(40,"CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT".length()), Read(2, 200, "8S8M17S"))))[0],
				onNegative(withReadName("0#0#34#readname", withQual(B("1"), withSequence("G", Read(0, 40, "1M")))))[0]
				));
		assertEquals(new BreakpointSummary(0, FWD, 5, 2, BWD, 200), e.getBreakendSummary());
		assertEquals("ATCGCAAGAGCG", e.getUntemplatedSequence());
		List<SAMRecordAssemblyEvidence> rl = e.getSubsequentRealignments();
		assertEquals(2, rl.size());
		
		RealignedSAMRecord r0 = (RealignedSAMRecordAssemblyEvidence)rl.get(0);
		assertEquals(new BreakpointSummary(2, FWD, 207, 1, BWD, 100), r0.getBreakendSummary());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTAT", S(r0.getAnchorSequence()));
		assertEquals("TCGAC", r0.getUntemplatedSequence());
		assertEquals(be.getEvidenceID() + "_0", r0.getEvidenceID());
		
		RealignedSAMRecord r1 = (RealignedSAMRecordAssemblyEvidence)rl.get(1);
		assertEquals(new BreakpointSummary(1, FWD, 104, 0, FWD, 40), r1.getBreakendSummary());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAA", S(r1.getAnchorSequence()));
		assertEquals("GTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", S(r1.getBreakendSequence()));
		assertEquals("GTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", SequenceUtil.reverseComplement(S(r1.getRemoteSAMRecord().getReadBases())));
		assertEquals("GTCA", r1.getUntemplatedSequence());
		assertEquals(be.getEvidenceID() + "_1", r1.getEvidenceID());
	}
	@Test
	public void getAllRealignments_should_return_all_breakpoints_bwd4() {
		GenomicProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES();
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(pc, aes, BWD, null, 0, 70, 10, B("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"), B("00000000000000000000000000000000000000000000000000"));
		//          1         2         3         4         5         6         7         8         9         0
		// 123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
		//          MMMMMMMMMM          MMMMMMMMMM          MMMMMMMMMM          MMMMMMMMMM          MMMMMMMMMM
		//                                                  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMM
		//                    SSSSSSSSSSMMMMMMMMMMSSSSSSSSSSSSSSSSSSSS offset = 0
		//          MMMMMMMMMM offset = 0
		//                                                            SSSSSSSSSSMMMMMMMMMM offset = 20
		//                                                  MMMMMMMMMM offset = 20                  
		RealignedSAMRecord re = (RealignedSAMRecordAssemblyEvidence)AssemblyFactory.incorporateRealignment(getContext(), e, ImmutableList.of(
				withSequence("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", withReadName("0#70#0#ReadName", Read(0, 30, "10S10M20S")))[0],
				withSequence("NNNNNNNNNN", withReadName("0#70#0#ReadName", Read(0, 10, "10M")))[0],
				withSequence("NNNNNNNNNNNNNNNNNNNN", withReadName("0#70#20#ReadName", Read(0, 70, "10S10M")))[0],
				withSequence("NNNNNNNNNN", withReadName("0#70#20#ReadName", Read(0, 50, "10M")))[0]
				));
		assertEquals(new BreakpointSummary(0, BWD, 70, 0, FWD, 79), re.getBreakendSummary());
		List<SAMRecordAssemblyEvidence> list = re.getSubsequentRealignments();
		Collections.sort(list, DirectedEvidenceOrder.ByStartEnd);
		assertEquals(new BreakpointSummary(0, BWD, 30, 0, FWD, 19), list.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(0, BWD, 50, 0, FWD, 39), list.get(1).getBreakendSummary());
		assertEquals(new BreakpointSummary(0, BWD, 70, 0, FWD, 59), list.get(2).getBreakendSummary());
	}
	*/
	@Test
	public void realign_should_shift_breakend_to_match_reference() {
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				0, 5, 5, B("AAAATTTT"), new byte[] {1,2,3,4,1,2,3,4});
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("TTTT", S(asEvidence(e).getBreakendSequence()));
		assertEquals(new BreakendSummary(0, FWD, 4), asEvidence(e).getBreakendSummary());
		assertEquals("4M4S", e.getCigarString());
	}
	@Test
	public void realign_should_align_to_reference_with_50bp_margin_around_expected_anchor_interval() {
		int margin = 50;
		for (int startpos = 300 - margin; startpos <= 300 + margin; startpos++) {
			String seq = S("N", 50) + S(Arrays.copyOfRange(RANDOM, 299, 399)); // genomic positions 300-400
			SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null,
					2, startpos, 100, B(seq), B(40, seq.length()));
			e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
			assertEquals(300, asEvidence(e).getBreakendSummary().start);
			assertEquals(50, asEvidence(e).getBreakendSequence().length);
		}
		// FWD breakend
		for (int startpos = 300 - margin; startpos <= 300 + margin; startpos++) {
			String seq = S(Arrays.copyOfRange(RANDOM, 299, 399)) + S("N", 50);
			SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
					2, startpos + 100 - 1, 100, B(seq), B(40, seq.length()));
			e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
			assertEquals(399, asEvidence(e).getBreakendSummary().start);
			assertEquals(50, asEvidence(e).getBreakendSequence().length);
		}
	}
	@Test
	public void realign_should_expand_window_by_breakend_length_to_allow_for_mapping_over_small_indels() {
		int indelSize = 20;
		String seq = "N" + S(Arrays.copyOfRange(RANDOM, 299-indelSize-100, 299-indelSize)) + S(Arrays.copyOfRange(RANDOM, 299, 399)); // genomic positions 300-400
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null,
				2, 300, 100, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("1S100M20D100M", e.getCigarString());
		
		seq = S(Arrays.copyOfRange(RANDOM, 299, 399)) + S(Arrays.copyOfRange(RANDOM, 399+indelSize, 399+indelSize+100)) + "N";
		e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				2, 399, 100, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("100M20D100M1S", e.getCigarString());
	}
	@Test
	public void realign_should_allow_small_anchor_deletion() {
		String seq = S(B('N', 100)) + S(Arrays.copyOfRange(RANDOM, 0, 100)) + S(Arrays.copyOfRange(RANDOM, 110, 210));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null,
				2, 1, 210, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("100S100M10D100M", e.getCigarString());
	}
	@Test
	public void realign_should_allow_small_anchor_insertion() {
		String seq = S(B('N', 100)) + S(Arrays.copyOfRange(RANDOM, 0, 100)) + "NNNNNNNNNN" + S(Arrays.copyOfRange(RANDOM, 100, 200));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null,
				2, 1, 200, B(seq), B(40, seq.length()));
		e = SAMRecordUtil.realign(getContext().getReference(), e, 50, true);
		assertEquals("100S100M10I100M", e.getCigarString());
	}
	@Test
	@Ignore("Not part of SAMRecordUtil realign functionality")
	public void realign_should_abort_if_anchor_turns_into_soft_clip() {
		String seq = S(Arrays.copyOfRange(RANDOM, 0, 10)) + S(Arrays.copyOfRange(RANDOM, 30, 70));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, null,
				2, 1, 10, B(seq), B(40, seq.length()));
		assertEquals("10M40S", e.getCigarString());
		assertEquals("10M40S", SAMRecordUtil.realign(getContext().getReference(), e, 50, true).getCigarString());
	}
	@Test
	@Ignore("Not part of SAMRecordUtil realign functionality")
	public void realign_should_not_flip_tandem_duplication() {
		ProcessingContext pc = getContext();
		String seq = S(Arrays.copyOfRange(RANDOM, 150, 200)) + S(Arrays.copyOfRange(RANDOM, 100, 200)) + S(Arrays.copyOfRange(RANDOM, 100, 150));
		SAMRecord e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), BWD, null,
				2, 100, 50, B(seq), B(40, seq.length()));
		
		// breakend alignment is actually better
		SAMRecord e1 = SAMRecordUtil.realign(getContext().getReference(), e, 200, true);
		assertEquals("48S102M50S", e1.getCigarString()); // 50S100M50S + homology
		
		// but if we constrain to require anchor bases, we shouldn't realign
		//SAMRecord e2 = SAMRecordUtil.realign(getContext().getReference(), e, 200, true, 0.5f); // 0.5f;
		// assertEquals("150S50M", e2.getCigarString());
	}
	@Test
	public void realign_should_turn_reference_bubble_into_reference_assembly() {
		SAMRecord ass = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), null,
				0, 10, 1,
				0, 17, 1,
				B("AAAAAAAA"),
				B("AAAAAAAA"));
		ass = SAMRecordUtil.realign(getContext().getReference(), ass, 50, true);
		assertEquals(0, SingleReadEvidence.createEvidence(SES(), ass).size());
	}
}

