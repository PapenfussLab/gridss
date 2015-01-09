package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.Sets;


public class SAMRecordAssemblyEvidenceTest extends TestHelper {
	@Test
	public void should_create_SAMRecord_for_assembly() {
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2);
		assertNotNull(e.getSAMRecord());
		assertEquals("GTAC", S(e.getSAMRecord().getReadBases()));
		assertArrayEquals( new byte[] {1,2,3,4}, e.getSAMRecord().getBaseQualities());
	}
	@Test
	public void should_create_placeholder_paired_read() {
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2);
		assertPairing(e.getSAMRecord(), e.getRemoteSAMRecord());
		assertTrue(e.getRemoteSAMRecord().getReadUnmappedFlag());
	}
	@Test
	public void anchor_positions_should_match_genomic() {
		SAMRecordAssemblyEvidence fwd = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(), 1, 10, 3, B("GTACCCA"), new byte[] { 1, 2, 3, 4, 4, 8 }, 2, 0);
		SAMRecordAssemblyEvidence bwd = AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(), 1, 10, 3, B("GTACCCA"), new byte[] { 1, 2, 3, 4, 4, 8 }, 2, 0);
		assertEquals(1, (int)fwd.getSAMRecord().getReferenceIndex());
		assertEquals(1, (int)bwd.getSAMRecord().getReferenceIndex());
		assertEquals(8, fwd.getSAMRecord().getAlignmentStart());
		assertEquals(10, bwd.getSAMRecord().getAlignmentStart());
		assertEquals("3M4S", fwd.getSAMRecord().getCigarString());
		assertEquals("4S3M", bwd.getSAMRecord().getCigarString());
	}
	@Test
	public void unanchor_positions_should_match_genomic() {
		SAMRecordAssemblyEvidence fwd = AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(1, FWD, 5, 10))), B("AAA"), B("AAA"), 2, 0);
		SAMRecordAssemblyEvidence bwd = AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(1, BWD, 5, 10))), B("AAA"), B("AAA"), 2, 0);
		assertEquals(1, (int)fwd.getSAMRecord().getReferenceIndex());
		assertEquals(1, (int)bwd.getSAMRecord().getReferenceIndex());
		assertEquals(5, fwd.getSAMRecord().getAlignmentStart());
		assertEquals(5, bwd.getSAMRecord().getAlignmentStart());
		assertEquals(10, fwd.getSAMRecord().getAlignmentEnd());
		assertEquals(10, bwd.getSAMRecord().getAlignmentEnd());
		assertEquals("1X4N1X3S", fwd.getSAMRecord().getCigarString());
		assertEquals("3S1X4N1X", bwd.getSAMRecord().getCigarString());
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
		assertEquals("1M3S", AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2).getSAMRecord().getCigarString());
		assertEquals("3S1M", AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
				1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2).getSAMRecord().getCigarString());
	}
	@Test
	public void should_use_XNXS_for_unanchored_interval_breakend() {
		assertEquals("1X297N1X4S", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
				NRRP(OEA(0, 1000, "1M", true))), B("GTAC"), new byte[] {1,2,3,4}, 0, 0).getSAMRecord().getCigarString());
		assertEquals("4S1X297N1X", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
				NRRP(OEA(0, 1000, "1M", false))), B("GTAC"), new byte[] {1,2,3,4}, 0, 0).getSAMRecord().getCigarString());
	}
	@Test
	public void should_use_minimal_cigar_representation() {
		assertEquals("1X3S", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(1, FWD, 5, 5))), B("AAA"), B("AAA"), 2, 0).getSAMRecord().getCigarString());
		assertEquals("2X3S", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(1, FWD, 5, 6))), B("AAA"), B("AAA"), 2, 0).getSAMRecord().getCigarString());
		assertEquals("1X1N1X3S", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(new MockDirectedEvidence(new BreakendSummary(1, FWD, 5, 7))), B("AAA"), B("AAA"), 2, 0).getSAMRecord().getCigarString());
	}
	@Test
	public void breakend_round_trip_should_be_unchanged() {
		for (SAMRecordAssemblyEvidence e : new SAMRecordAssemblyEvidence[] {
				AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(NRRP(OEA(0, 1000, "1M", true))), B("GTAC"), new byte[] {1,2,3,4}, 0, 0),
				AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(NRRP(OEA(0, 1000, "1M", false))), B("GTAC"), new byte[] {1,2,3,4}, 0, 0),
				AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(), 1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2),
				AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(), 1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2),
				big(),
			}) {
			SAMRecordAssemblyEvidence r = new SAMRecordAssemblyEvidence(e.getEvidenceSource(), e.getSAMRecord(), null);
			assertEvidenceEquals(e, r);
		}
	}
	public void assertEvidenceEquals(AssemblyEvidence e, AssemblyEvidence r) {
		assertEquals(e.getAssemblyAnchorLength(), r.getAssemblyAnchorLength());
		//assertArrayEquals(e.getAssemblyAnchorQuals() , r.getAssemblyAnchorQuals());
		assertEquals(S(e.getAssemblyAnchorSequence()) , S(r.getAssemblyAnchorSequence()));
		assertEquals(e.getAssemblyBaseCount(EvidenceSubset.NORMAL) , r.getAssemblyBaseCount(EvidenceSubset.NORMAL));
		assertEquals(e.getAssemblyBaseCount(EvidenceSubset.TUMOUR) , r.getAssemblyBaseCount(EvidenceSubset.TUMOUR));
		assertEquals(e.getAssemblyBaseCount(EvidenceSubset.ALL) , r.getAssemblyBaseCount(EvidenceSubset.ALL));
		assertEquals(e.getAssemblyReadPairLengthMax(EvidenceSubset.NORMAL) , r.getAssemblyReadPairLengthMax(EvidenceSubset.NORMAL));
		assertEquals(e.getAssemblyReadPairLengthMax(EvidenceSubset.TUMOUR) , r.getAssemblyReadPairLengthMax(EvidenceSubset.TUMOUR));
		assertEquals(e.getAssemblyReadPairLengthMax(EvidenceSubset.ALL) , r.getAssemblyReadPairLengthMax(EvidenceSubset.ALL));
		assertEquals(S(e.getAssemblySequence()) , S(r.getAssemblySequence()));
		assertEquals(e.getAssemblySoftClipLengthMax(EvidenceSubset.NORMAL) , r.getAssemblySoftClipLengthMax(EvidenceSubset.NORMAL));
		assertEquals(e.getAssemblySoftClipLengthMax(EvidenceSubset.TUMOUR) , r.getAssemblySoftClipLengthMax(EvidenceSubset.TUMOUR));
		assertEquals(e.getAssemblySoftClipLengthMax(EvidenceSubset.ALL) , r.getAssemblySoftClipLengthMax(EvidenceSubset.ALL));
		assertEquals(e.getAssemblySoftClipLengthTotal(EvidenceSubset.NORMAL) , r.getAssemblySoftClipLengthTotal(EvidenceSubset.NORMAL));
		assertEquals(e.getAssemblySoftClipLengthTotal(EvidenceSubset.TUMOUR) , r.getAssemblySoftClipLengthTotal(EvidenceSubset.TUMOUR));
		assertEquals(e.getAssemblySoftClipLengthTotal(EvidenceSubset.ALL) , r.getAssemblySoftClipLengthTotal(EvidenceSubset.ALL));
		assertEquals(e.getAssemblySupportCountReadPair(EvidenceSubset.NORMAL) , r.getAssemblySupportCountReadPair(EvidenceSubset.NORMAL));
		assertEquals(e.getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR) , r.getAssemblySupportCountReadPair(EvidenceSubset.TUMOUR));
		assertEquals(e.getAssemblySupportCountReadPair(EvidenceSubset.ALL) , r.getAssemblySupportCountReadPair(EvidenceSubset.ALL));
		assertEquals(e.getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL) , r.getAssemblySupportCountSoftClip(EvidenceSubset.NORMAL));
		assertEquals(e.getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR) , r.getAssemblySupportCountSoftClip(EvidenceSubset.TUMOUR));
		assertEquals(e.getAssemblySupportCountSoftClip(EvidenceSubset.ALL) , r.getAssemblySupportCountSoftClip(EvidenceSubset.ALL));
		assertArrayEquals(e.getBreakendQuality() , r.getBreakendQuality());
		assertEquals(S(e.getBreakendSequence()) , S(r.getBreakendSequence()));
		assertEquals(e.getBreakendSummary(), r.getBreakendSummary());
		assertEquals(e.getBreakendSummary().getClass(), r.getBreakendSummary().getClass());
		assertEquals(e.getEvidenceID() , r.getEvidenceID());
		assertEquals(e.getEvidenceSource() , r.getEvidenceSource());
		assertEquals(e.getFilters() , r.getFilters());
		assertEquals(e.getLocalBaseLength() , r.getLocalBaseLength());
		assertEquals(e.getLocalMapq() , r.getLocalMapq());
		assertEquals(e.getLocalMaxBaseQual() , r.getLocalMaxBaseQual());
		assertEquals(e.getLocalTotalBaseQual() , r.getLocalTotalBaseQual());
		if (e instanceof DirectedBreakpoint) {
			assertTrue(r instanceof DirectedBreakpoint);
			DirectedBreakpoint de = (DirectedBreakpoint)e;
			DirectedBreakpoint dr = (DirectedBreakpoint)r;
			assertEquals(de.getBreakendSummary(), dr.getBreakendSummary());
			assertEquals(de.getRemoteMapq(), dr.getRemoteMapq());
			assertEquals(de.getRemoteBaseLength(), dr.getRemoteBaseLength());
			assertEquals(de.getRemoteBaseCount(), dr.getRemoteBaseCount());
			assertEquals(de.getRemoteMaxBaseQual(), dr.getRemoteMaxBaseQual());
			assertEquals(de.getRemoteTotalBaseQual(), dr.getRemoteTotalBaseQual());
			assertEquals(de.getUntemplatedSequence(), dr.getUntemplatedSequence());
		}
	}
	private SAMRecordAssemblyEvidence big() {
		return new AssemblyFactoryTest().big();
	}
	@Test
	public void should_track_breakend_evidence() {
		DirectedEvidence e1 = SCE(BWD, Read(0, 1, "5S5M"));
		DirectedEvidence e2 = SCE(BWD, Read(0, 1, "6S5M"));
		DirectedEvidence e3 = NRRP(OEA(0, 1, "1M", false));
		DirectedEvidence b1 = SCE(BWD, Read(0, 1, "7S5M"));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(e1, e2, e3),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2);
		assertTrue(e.isPartOfAssemblyBreakend(e1));
		assertTrue(e.isPartOfAssemblyBreakend(e2));
		assertTrue(e.isPartOfAssemblyBreakend(e3));
		assertFalse(e.isPartOfAssemblyBreakend(b1));
		e = new SAMRecordAssemblyEvidence(AES(), e.getSAMRecord(), null);
		assertTrue(e.isPartOfAssemblyBreakend(e1));
		assertTrue(e.isPartOfAssemblyBreakend(e2));
		assertTrue(e.isPartOfAssemblyBreakend(e3));
		assertFalse(e.isPartOfAssemblyBreakend(b1));
	}
	@Test
	public void should_rehydrate_breakend_evidence() {
		DirectedEvidence e1 = SCE(BWD, Read(0, 1, "5S5M"));
		DirectedEvidence e2 = SCE(BWD, Read(0, 1, "6S5M"));
		DirectedEvidence e3 = NRRP(OEA(0, 1, "1M", false));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(e1, e2, e3),
			1, 2, 1, B("GTAC"), new byte[] {1,2,3,4}, 1, 2);
		e.hydrateEvidenceSet(e1);
		e.hydrateEvidenceSet(e2);
		e.hydrateEvidenceSet(e3);
		assertEquals(3, e.getEvidence().size());
		e = new SAMRecordAssemblyEvidence(AES(), e.getSAMRecord(), null);
		e.hydrateEvidenceSet(e1);
		e.hydrateEvidenceSet(e2);
		e.hydrateEvidenceSet(e3);
		assertEquals(3, e.getEvidence().size());
	}
	@Test
	public void realign_should_shift_breakend_to_match_reference() {
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 5, 5, B("AAAATTTT"), new byte[] {1,2,3,4,1,2,3,4}, 0, 0).realign();
		assertEquals("TTTT", S(e.getBreakendSequence()));
		assertEquals(new BreakendSummary(0, FWD, 4, 4), e.getBreakendSummary());
		assertEquals("4M4S", e.getSAMRecord().getCigarString());
	}
	@Test
	public void realign_should_align_to_reference_with_30bp_margin() {
		int margin = 30;
		for (int startpos = 41 - margin; startpos <= 41 + margin; startpos++) {
			SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
					2, startpos, 80-40, B(S(B('N', 60)) + S(Arrays.copyOfRange(RANDOM, 40, 80))), B(40, 100), 0, 0).realign();
			assertEquals(41, e.getBreakendSummary().start);
			assertEquals(60, e.getBreakendSequence().length);
		}
	}
	@Test
	public void realign_should_allow_small_anchor_deletion() {
		String seq = S(B('N', 100)) + S(Arrays.copyOfRange(RANDOM, 0, 100)) + S(Arrays.copyOfRange(RANDOM, 110, 210));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
				2, 1, 210, B(seq), B(40, seq.length()), 0, 0).realign();
		assertEquals("100S100M10D100M", e.getSAMRecord().getCigarString());
	}
	@Test
	public void realign_should_allow_small_anchor_insertion() {
		String seq = S(B('N', 100)) + S(Arrays.copyOfRange(RANDOM, 0, 100)) + "NNNNNNNNNN" + S(Arrays.copyOfRange(RANDOM, 100, 200));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
				2, 1, 200, B(seq), B(40, seq.length()), 0, 0).realign();
		assertEquals("100S100M10I100M", e.getSAMRecord().getCigarString());
	}
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	@Test
	public void realign_should_fix_778_chr1_170849702_misalignment() {
		String assembly = "ATCCATCCCTATGACCCAAACATCTCCCACCAGGCCTCATGTTCAATATTAAAGATCACATTTCAACTTGAGATTTGGAGGGGACAAACATACAAATCATATCATTATCTCTCTCCCCACTTCTCTCTTTATCAATCCCTCCCTCTTTGTCAATCTTAGCCTTGGCCTTCAGATTTTACCACTTGATTTTTCACATTTTCTGTATTCTTAAT"
				+ "GATTATTATATTTTCATGTTCTTGCTAATCTATATCATGGTTAGAAATCAAAGCATGCCGAAATTTCTCTCTTACTTTTTTTGCTGTT";
		File ref = new File("src/test/resources/chr1_170849600_170849850.fa");
		ProcessingContext context = new ProcessingContext(getFSContext(),
				new ArrayList<Header>(), new SoftClipParameters(),
				new ReadPairParameters(), new AssemblyParameters(),
				new RealignmentParameters(), new VariantCallingParameters(),
				ref, false, false);
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(context, AES(context), BWD, Sets.<DirectedEvidence>newHashSet(),
				0, 170849702-170849600+1, 97, B(assembly), B(40, assembly.length()), 0, 0).realign();
		// anchor location is 11bp off
		assertEquals("212S88M", e.getSAMRecord().getCigarString());
		assertEquals(170849713-170849600+1, e.getBreakendSummary().start);
	}
	@Test
	public void realign_should_be_abandoned_if_no_soft_clip_left() {
		String assembly = "AAAAAAAAAAA";
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 5, B(assembly), B(40, assembly.length()), 0, 0);
		SAMRecordAssemblyEvidence re = e.realign();
		assertEquals(e, re);
	}
}