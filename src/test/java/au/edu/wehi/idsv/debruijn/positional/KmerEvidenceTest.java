package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.stream.IntStream;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;


public class KmerEvidenceTest extends TestHelper {
	@Test
	public void softclip() {
		for (SoftClipEvidence sce : new SoftClipEvidence[] {
				SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1S4M6S")))),
				SCE(BWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1S4M6S"))))}) {
			KmerEvidence e = KmerEvidence.create(4, sce, false);
			assertEquals(1, e.startPosition());
			assertEquals(1, e.endPosition());
			assertEquals(11-4+1, e.length());
			assertEquals(1, e.weight(0));
			assertEquals(2, e.weight(1));
			assertEquals(3, e.weight(2));
			assertEquals(4, e.weight(3));
			assertEquals(5, e.weight(4));
			assertEquals(6, e.weight(5));
			assertEquals(7, e.weight(6));
			assertEquals(8, e.weight(7));			
			assertEquals("ACGT", K(4, e.kmer(0)));
			assertEquals("CGTT", K(4, e.kmer(1)));
			assertEquals("GTTA", K(4, e.kmer(2)));
			assertEquals("TTAT", K(4, e.kmer(3)));
			assertEquals("TATA", K(4, e.kmer(4)));
			assertEquals("ATAC", K(4, e.kmer(5)));
			assertEquals("TACC", K(4, e.kmer(6)));
			assertEquals("ACCG", K(4, e.kmer(7)));
			assertFalse(e.isAnchored(0));
			assertTrue(e.isAnchored(1));
			assertFalse(e.isAnchored(2));
			assertFalse(e.isAnchored(3));
			assertFalse(e.isAnchored(4));
			assertFalse(e.isAnchored(5));
			assertFalse(e.isAnchored(6));
			assertFalse(e.isAnchored(7));
			assertEquals(sce.getEvidenceID(), e.evidence().getEvidenceID());
			for (int i = 0; i < e.length(); i++) {
				assertEquals(1 + i, e.node(i).lastStart());
				assertEquals(1 + i, e.node(i).lastEnd());
				assertEquals(e.isAnchored(i), e.node(i).isReference());
				assertEquals(e.kmer(i), e.node(i).lastKmer());
				assertEquals(e.weight(i), e.node(i).weight());
				assertEquals(e, e.node(i).evidence());
			}
		}
	}
	@Test
	public void pair_fwd() {
		//          1         2         3
		// 123456789012345678901234567890123456789012345678901234567890
		// SMMMMMMMMMS                
		// |------------------| min
		//          MMMMMMMMMMM
		// |----------------------------| max
		//                    MMMMMMMMMMM
		MockSAMEvidenceSource ses = SES(20, 30);
		KmerEvidence e = KmerEvidence.create(4, NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 2, "1S9M1S", true, 1, 1, "11M", false)))));
		assertEquals(10, e.startPosition());
		assertEquals(20, e.endPosition());
		assertEquals(11-4+1, e.length());
		assertEquals("ACGT", K(4, e.kmer(0)));
		assertEquals("CGTT", K(4, e.kmer(1)));
		assertEquals("GTTA", K(4, e.kmer(2)));
		assertEquals("TTAT", K(4, e.kmer(3)));
		assertEquals("TATA", K(4, e.kmer(4)));
		assertEquals("ATAC", K(4, e.kmer(5)));
		assertEquals("TACC", K(4, e.kmer(6)));
		assertEquals("ACCG", K(4, e.kmer(7)));
		assertEquals(1, e.weight(0));
		assertEquals(2, e.weight(1));
		assertEquals(3, e.weight(2));
		assertEquals(4, e.weight(3));
		assertEquals(5, e.weight(4));
		assertEquals(6, e.weight(5));
		assertEquals(7, e.weight(6));
		assertEquals(8, e.weight(7));
		for (int i = 0; i < e.length(); i++) assertFalse(e.isAnchored(i));
	}
	@Test
	public void pair_bwd() {
		//          1         2         3
		// 123456789012345678901234567890123456789012345678901234567890
		//                    SMMMMMMMMMS                
		//           |------------------| min
		//           MMMMMMMMMMM
		// |----------------------------| max
		// MMMMMMMMMMM
		MockSAMEvidenceSource ses = SES(20, 30);
		KmerEvidence e = KmerEvidence.create(4, NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", true)))));
		assertEquals(1, e.startPosition());
		assertEquals(11, e.endPosition());
		assertEquals(11-4+1, e.length());
		assertEquals("ACGT", K(4, e.kmer(0)));
		assertEquals("CGTT", K(4, e.kmer(1)));
		assertEquals("GTTA", K(4, e.kmer(2)));
		assertEquals("TTAT", K(4, e.kmer(3)));
		assertEquals("TATA", K(4, e.kmer(4)));
		assertEquals("ATAC", K(4, e.kmer(5)));
		assertEquals("TACC", K(4, e.kmer(6)));
		assertEquals("ACCG", K(4, e.kmer(7)));
		assertEquals(1, e.weight(0));
		assertEquals(2, e.weight(1));
		assertEquals(3, e.weight(2));
		assertEquals(4, e.weight(3));
		assertEquals(5, e.weight(4));
		assertEquals(6, e.weight(5));
		assertEquals(7, e.weight(6));
		assertEquals(8, e.weight(7));
		for (int i = 0; i < e.length(); i++) assertFalse(e.isAnchored(i));
	}
	@Test
	public void should_reverse_comp_pair_if_required() {
		MockSAMEvidenceSource ses = SES(20, 30);
		KmerEvidence e = KmerEvidence.create(4, NRRP(ses, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
		assertEquals("CGGT", K(4, e.kmer(0))); // reverse comp of ending ACCG
	}
	@Test
	public void softclip_should_trim_soft_clip_on_other_side() {
		assertEquals(4, KmerEvidence.create(2, SCE(FWD, Read(0, 1, "10H1S2M3S")), true).length());
		assertEquals(2, KmerEvidence.create(2, SCE(BWD, Read(0, 1, "10H1S2M3S")), true).length());
		assertEquals(7, KmerEvidence.create(2, SCE(FWD, Read(0, 1, "5M3S")), true).length());
		assertEquals(7, KmerEvidence.create(2, SCE(BWD, Read(0, 1, "3S5M")), true).length());
		
	}
	@Test
	public void should_return_null_if_kmer_longer_than_read() {
		assertNotNull(KmerEvidence.create(12, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S"))), true));
		assertNull(KmerEvidence.create(13, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S"))), true));
		assertNotNull(KmerEvidence.create(11, NRRP(withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
		assertNull(KmerEvidence.create(12, NRRP(withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
	}
	@Test
	public void short_sc_anchor_should_be_considered_unanchored_evidence() {
		KmerEvidence e = KmerEvidence.create(3, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "2M10S"))), true);
		assertTrue(IntStream.range(0, e.length()).allMatch(i -> !e.isAnchored(i)));
	}
	@Test
	public void should_exclude_ambiguous_kmers() {
		KmerEvidence e = KmerEvidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTNATACCG", Read(0, 2, "1S4M6S")))), false);
		assertNotNull(e.node(0));
		assertNull(e.node(1));
		assertNull(e.node(2));
		assertNull(e.node(3));
		assertNull(e.node(4));
		assertNotNull(e.node(5));
		e = KmerEvidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4}, withSequence("ACGTN", Read(0, 2, "4M1S")))), false);
		assertNotNull(e.node(0));
		assertNull(e.node(1));
		e = KmerEvidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("NCGTTATACCN", Read(0, 2, "1S4M6S")))), false);
		assertNull(e.node(0));
		assertNotNull(e.node(1));
	}
	@Test
	@Ignore
	public void should_handle_cigar_indels() {
		// TODO: make sure both sides are at the correct position for cigars such as 5S1M1D1M5S
	}
	@Test
	public void createAnchor_should_anchor_at_alignment_closest_to_breakend() {
		assertEquals(100, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1M2I3M", false)), 0, null).startPosition());
		assertEquals(99, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1S1M2I3M", false)), 0, null).startPosition());
		assertEquals(98, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "2S1M2I3M", false)), 0, null).startPosition());
		assertEquals(98, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1M2I3M", true)), 0, null).startPosition());
		assertEquals(98, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1M2I3M1S", true)), 0, null).startPosition());
		assertEquals(98, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1M2I3M2S", true)), 0, null).startPosition());
		assertEquals(97, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1S1M2I3M", true)), 0, null).startPosition());
		assertEquals(97, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "1S1M2I3M2S", true)), 0, null).startPosition());
		assertEquals(96, KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "2S1M2I3M10S", true)), 0, null).startPosition());
	}
	@Test
	public void createAnchor_should_support_anchor() {
		KmerEvidence a = KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "10M2I30M", false)), 0, null);
		for (int i = 0; i < a.length(); i++) {
			if (a.node(i) != null) {
				assertTrue(a.node(i).isReference());
			}
		}
	}
	@Test
	public void createAnchor_should_only_anchor_subalignment_closest_to_breakend() {
		KmerEvidence a = KmerEvidence.createAnchor(1, NRRP(SES(), OEA(0, 100, "10M2I4M", true)), 0, null);
		for (int i = 0; i < 12; i++) {
			assertNull(a.node(i));
		}
		a = KmerEvidence.createAnchor(4, NRRP(SES(), OEA(0, 100, "10M2I4M", false)), 0, null);
		for (int i = 10-3+1; i < a.length(); i++) {
			assertNull(a.node(i));
		}
	}
	private int nodesNotNull(KmerEvidence a) {
		int n = 0;
		for (int i = 0; i < a.length(); i++) {
			if (a.node(i) != null) n++;
		}
		return n;
	}
	@Test
	public void createAnchor_should_not_include_mismatches_near_end_of_read() {
		assertEquals(10, nodesNotNull(KmerEvidence.createAnchor(1, NRRP(SES(), withSequence("TTAAAAAAAA", OEA(0, 100, "10M", true))), 0, SMALL_FA)));
		assertEquals(9, nodesNotNull(KmerEvidence.createAnchor(1, NRRP(SES(), withSequence("TTAAAAAAAA", OEA(0, 100, "10M", true))), 1, SMALL_FA)));
		assertEquals(8, nodesNotNull(KmerEvidence.createAnchor(1, NRRP(SES(), withSequence("TTAAAAAAAA", OEA(0, 100, "10M", true))), 2, SMALL_FA)));
		assertEquals(8, nodesNotNull(KmerEvidence.createAnchor(1, NRRP(SES(), withSequence("TTAAAAAAAA", OEA(0, 100, "10M", true))), 3, SMALL_FA)));
	}
	@Test
	public void createAnchor_should_check_reference_overrun() {
		SoftClipEvidence placholderEvidence = SCE(BWD, Read(0, 1, "2S5M2S"));
		assertEquals(5, nodesNotNull(KmerEvidence.createAnchor(placholderEvidence, 1, Read(0, 1, "2S5M2S"), BWD, 2, SMALL_FA)));
		assertEquals(5, nodesNotNull(KmerEvidence.createAnchor(placholderEvidence, 1, Read(0, 1, "2S5M2S"), FWD, 2, SMALL_FA)));
		// overruns on both end
		assertEquals(2, nodesNotNull(KmerEvidence.createAnchor(placholderEvidence, 1, Read(0, -10, "10M"), FWD, 4, SMALL_FA)));
		// end overrun
		assertEquals(1, nodesNotNull(KmerEvidence.createAnchor(placholderEvidence, 1, Read(0, 10000, "10M"), FWD, 10, SMALL_FA)));
	}
}
