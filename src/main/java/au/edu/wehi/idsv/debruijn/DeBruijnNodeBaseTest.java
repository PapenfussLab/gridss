package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.*;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.TestHelper;


public class DeBruijnNodeBaseTest extends TestHelper {
	private static VariantEvidence fsc = new VariantEvidence(4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))), getContext().getLinear());
	private static VariantEvidence bsc = new VariantEvidence(4, SCE(BWD, withSequence("TTAACCGGCCAATT", Read(0, 17, "7S7M"))), getContext().getLinear());
	private static VariantEvidence frp = new VariantEvidence(4, NRRP(DP(0, 1, "10M", true, 1, 1, "10M", false)), getContext().getLinear());
	@Test
	public void getExpectedPosition_should_match_variant_evidence_expected_position() {
		for (int i = 0; i < 4; i++) {
			assertEquals(fsc.getExpectedLinearPosition(i), new DeBruijnNodeBase(fsc, i, new ReadKmer(0, 2, false)).getExpectedPosition(), 0);
			assertEquals(fsc.getExpectedLinearPosition(i), new DeBruijnNodeBase(fsc, i, new ReadKmer(0, 2, false)).getExpectedPosition(), 0);
		}
	}
	@Test
	public void getExpectedPosition_should_be_weighted_mean_of_evidence_expected_position() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 1, new ReadKmer(0, 1, false));
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 4, new ReadKmer(0, 2, false));
		assertEquals(LCCB + 11, r1.getExpectedPosition(), 0);
		assertEquals(LCCB + 14, r2.getExpectedPosition(), 0);
		r1.add(r2);
		assertEquals(LCCB + 13, r1.getExpectedPosition(), 0);
	}
	@Test
	public void should_track_reference_kmers() {
		int offset = 0;
		for (ReadKmer rk : fsc.getReadKmers()) {
			assertEquals(offset < 4, new DeBruijnNodeBase(fsc, offset, rk).isReference());
			offset++;
		}
		offset = 0;
		for (ReadKmer rk : fsc.getReadKmers()) {
			assertEquals(offset >= 7, new DeBruijnNodeBase(bsc, offset, rk).isReference());
			offset++;
		}
		for (ReadKmer rk : fsc.getReadKmers()) {
			assertFalse(new DeBruijnNodeBase(frp, offset, rk).isReference()); 
		}
	}
	@Test
	public void should_track_support() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 0, new ReadKmer(0, 1, false));
		assertEquals(1, r1.getSupportingEvidenceList().size());
		assertEquals(1, r2.getSupportingEvidenceList().size());
		assertTrue(r1.getSupportingEvidenceList().contains(fsc.getDirectedEvidence()));
		assertTrue(r2.getSupportingEvidenceList().contains(bsc.getDirectedEvidence()));
		r1.add(r2);
		assertEquals(2, r1.getSupportingEvidenceList().size());
		assertTrue(r1.getSupportingEvidenceList().contains(fsc.getDirectedEvidence()));
		assertTrue(r1.getSupportingEvidenceList().contains(bsc.getDirectedEvidence()));
	}
	@Test
	public void should_track_weight() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 0, new ReadKmer(0, 2, false));
		assertEquals(1, r1.getWeight());
		assertEquals(2, r2.getWeight());
		r1.add(r2);
		assertEquals(3, r1.getWeight());
	}
	@Test
	public void reference_with_non_reference_should_be_reference() {
		DeBruijnNodeBase r = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase nr = new DeBruijnNodeBase(fsc, 10, new ReadKmer(0, 1, false));
		assertTrue(r.isReference());
		assertFalse(nr.isReference());
		r.add(nr);
		assertTrue(r.isReference());
		
		r = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		nr = new DeBruijnNodeBase(fsc, 10, new ReadKmer(0, 1, false));
		nr.add(r);
		assertTrue(nr.isReference());
	}
	@Test
	public void getExpectedPositionForDirectAnchor_FWD_should_return_expected_position_immediately_before_first_kmer() {
		DeBruijnNodeBase r0 = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 1, new ReadKmer(0, 1, false));
		assertEquals(r0.getExpectedPosition() - 1, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(r0)), 0);
		assertEquals(r0.getExpectedPosition() - 1, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(r0, r1)), 0);
		assertEquals(r0.getExpectedPosition() - 1 - 0.5, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(r0, r0)), 0);
	}
	@Test
	public void getExpectedPositionForDirectAnchor_BWD_should_return_expected_position_immediately_before_first_kmer() {
		DeBruijnNodeBase r0 = new DeBruijnNodeBase(bsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(bsc, 1, new ReadKmer(0, 1, false));
		assertEquals(r0.getExpectedPosition() + 1, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BWD, ImmutableList.of(r0)), 0);
		assertEquals(r1.getExpectedPosition() + 1, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BWD, ImmutableList.of(r0, r1)), 0);
		assertEquals(r0.getExpectedPosition() + 1 + 0.5, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BWD, ImmutableList.of(r0, r0)), 0);
	}
	@Test
	public void getExpectedPositionForDirectAnchor_should_weight_expected_positions_by_kmer_weight() {
		DeBruijnNodeBase w1 = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase w2 = new DeBruijnNodeBase(fsc, 4, new ReadKmer(0, 2, false));
		assertEquals(LCCB + 10, w1.getExpectedPosition(), 0);
		assertEquals(LCCB + 14, w2.getExpectedPosition(), 0);
		// Anchor, w1, w2
		// w1 expected at 10 -> anchor at 9
		// w2 expected at 14 -> anchor at 12
		assertEquals(LCCB + 11, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(w1, w2)), 0);
	}
	@Test
	public void getExpectedPositionForDirectAnchor_should_only_consider_direct_evidence_in_given_direction() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 0, new ReadKmer(0, 1, false));
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 0, new ReadKmer(0, 2, false));
		DeBruijnNodeBase r3 = new DeBruijnNodeBase(frp, 0, new ReadKmer(0, 2, false));
		assertEquals(r1.getExpectedPosition() - 1, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(r1)), 0);
		assertTrue(Double.isNaN(DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BWD, ImmutableList.of(r1))));
		assertTrue(Double.isNaN(DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(r2))));
		assertEquals(r2.getExpectedPosition() + 1, DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BWD, ImmutableList.of(r2)), 0);
		assertTrue(Double.isNaN(DeBruijnNodeBase.getExpectedPositionForDirectAnchor(FWD, ImmutableList.of(r3))));
		assertTrue(Double.isNaN(DeBruijnNodeBase.getExpectedPositionForDirectAnchor(BWD, ImmutableList.of(r3))));
	}
}
