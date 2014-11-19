package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

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
		assertEquals(10, bwd.getSAMRecord().getAlignmentStart());
		assertEquals("1X5N3S", fwd.getSAMRecord().getCigarString());
		assertEquals("3S5N1X", bwd.getSAMRecord().getCigarString());
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
	public void should_use_XNS_for_unanchored_interval_breakend() {
		assertEquals("1X298N4S", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
				NRRP(OEA(0, 1000, "1M", true))), B("GTAC"), new byte[] {1,2,3,4}, 0, 0).getSAMRecord().getCigarString());
		assertEquals("4S298N1X", AssemblyFactory.createUnanchored(getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
				NRRP(OEA(0, 1000, "1M", false))), B("GTAC"), new byte[] {1,2,3,4}, 0, 0).getSAMRecord().getCigarString());
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
}