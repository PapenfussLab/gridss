package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class DeBruijnNodeBaseTest extends TestHelper {
	private static VariantEvidence fsc = new VariantEvidence(4, SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))), getContext().getLinear());
	private static VariantEvidence bsc = new VariantEvidence(4, SCE(BWD, withSequence("TTAACCGGCCAATT", Read(0, 17, "7S7M"))), getContext().getLinear());
	private static VariantEvidence frp = new VariantEvidence(4, NRRP(DP(0, 1, "10M", true, 1, 1, "10M", false)), getContext().getLinear());
	@Test
	public void getExpectedPosition_should_match_variant_evidence_expected_position() {
		for (int i = 0; i < 4; i++) {
			assertEquals(fsc.getExpectedLinearPosition(i), new DeBruijnNodeBase(fsc, i).getExpectedPosition(), 0);
			assertEquals(fsc.getExpectedLinearPosition(i), new DeBruijnNodeBase(fsc, i).getExpectedPosition(), 0);
		}
	}
	@Test
	public void getExpectedPosition_should_be_weighted_mean_of_evidence_expected_position() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 1);
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 4);
		assertEquals(LCCB + 11, r1.getExpectedPosition(), 0);
		assertEquals(LCCB + 14, r2.getExpectedPosition(), 0);
		r1.add(r2);
		assertEquals(LCCB + 13, r1.getExpectedPosition(), 0);
	}
	@Test
	public void should_track_reference_kmers() {
		for (int i = 0; i < fsc.getKmers().length(); i++) {
			assertEquals(i < 4, new DeBruijnNodeBase(fsc, i).isReference());
		}
		for (int i = 0; i < bsc.getKmers().length(); i++) {
			assertEquals(i >= 7, new DeBruijnNodeBase(bsc, i).isReference());
		}
		for (int i = 0; i < frp.getKmers().length(); i++) {
			assertFalse(new DeBruijnNodeBase(frp, i).isReference());
		}
	}
	@Test
	public void should_track_support_evidenceIDs() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 0);
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 0);
		assertEquals(1, r1.getSupportingEvidenceList().size());
		assertEquals(1, r2.getSupportingEvidenceList().size());
		assertTrue(r1.getSupportingEvidenceList().contains(fsc.getEvidenceID()));
		assertTrue(r2.getSupportingEvidenceList().contains(bsc.getEvidenceID()));
		r1.add(r2);
		assertEquals(2, r1.getSupportingEvidenceList().size());
		assertTrue(r1.getSupportingEvidenceList().contains(fsc.getEvidenceID()));
		assertTrue(r1.getSupportingEvidenceList().contains(bsc.getEvidenceID()));
	}
	@Test
	public void should_track_weight() {
		DeBruijnNodeBase r1 = new DeBruijnNodeBase(fsc, 0);
		DeBruijnNodeBase r2 = new DeBruijnNodeBase(bsc, 0);
		assertEquals(31, r1.getWeight());
		assertEquals(31, r2.getWeight());
		r1.add(r2);
		assertEquals(62, r1.getWeight());
	}
	@Test
	public void reference_with_non_reference_should_be_reference() {
		DeBruijnNodeBase r = new DeBruijnNodeBase(fsc, 0);
		DeBruijnNodeBase nr = new DeBruijnNodeBase(fsc, 10);
		assertTrue(r.isReference());
		assertFalse(nr.isReference());
		r.add(nr);
		assertTrue(r.isReference());
		
		r = new DeBruijnNodeBase(fsc, 0);
		nr = new DeBruijnNodeBase(fsc, 10);
		nr.add(r);
		assertTrue(nr.isReference());
	}
	@Test
	public void getExpectedPosition_should_return_max_weight_position() {
		VariantEvidence v = new VariantEvidence(1, SCE(FWD, withSequence("TTTTTTT", Read(0, 1, "1M6S"))), getContext().getLinear());
		DeBruijnNodeBase r = new DeBruijnNodeBase(v, 0);
		r.add(new DeBruijnNodeBase(v, 1));
		r.add(new DeBruijnNodeBase(v, 2));
		r.add(new DeBruijnNodeBase(v, 3));
		r.add(new DeBruijnNodeBase(v, 3));
		r.add(new DeBruijnNodeBase(v, 4));
		r.add(new DeBruijnNodeBase(v, 5));
		r.add(new DeBruijnNodeBase(v, 6));
		assertEquals(LCCB + 4, r.getExpectedPosition()); // position 3 has max total weight
	}
	@Test
	public void getExpectedReferencePosition_should_return_max_weight_reference_position() {
		VariantEvidence v = new VariantEvidence(1, SCE(FWD, withSequence("TTTT", Read(0, 1, "2M2S"))), getContext().getLinear());
		DeBruijnNodeBase r = new DeBruijnNodeBase(v, 0);
		r.add(new DeBruijnNodeBase(v, 0));
		r.add(new DeBruijnNodeBase(v, 0));
		r.add(new DeBruijnNodeBase(v, 1));
		r.add(new DeBruijnNodeBase(v, 2));
		r.add(new DeBruijnNodeBase(v, 3));
		assertEquals(LCCB + 1, r.getExpectedReferencePosition()); // 4 & 8 positions are non-reference
	}
}
