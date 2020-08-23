package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.*;
import org.apache.commons.lang3.NotImplementedException;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.stream.IntStream;

import static org.junit.Assert.*;


public class KmerEvidenceTest extends TestHelper {
	@Test
	public void fwd_softclip() {
		SoftClipEvidence sce = SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1S4M6S"))));
		KmerEvidence e = KmerEvidence.create(4, sce);
		//int startsclength = 1;
		assertEquals(2, e.startPosition());
		assertEquals(2, e.endPosition());
		assertEquals(11-1-(4-1), e.length());
		assertEquals(2, e.weight(0));
		assertEquals(3, e.weight(1));
		assertEquals(4, e.weight(2));
		assertEquals(5, e.weight(3));
		assertEquals(6, e.weight(4));
		assertEquals(7, e.weight(5));
		assertEquals(8, e.weight(6));	
		assertEquals("CGTT", K(4, e.kmer(0)));
		assertEquals("GTTA", K(4, e.kmer(1)));
		assertEquals("TTAT", K(4, e.kmer(2)));
		assertEquals("TATA", K(4, e.kmer(3)));
		assertEquals("ATAC", K(4, e.kmer(4)));
		assertEquals("TACC", K(4, e.kmer(5)));
		assertEquals("ACCG", K(4, e.kmer(6)));
		assertTrue(e.isAnchored(0));
		assertFalse(e.isAnchored(1));
		assertFalse(e.isAnchored(2));
		assertFalse(e.isAnchored(3));
		assertFalse(e.isAnchored(4));
		assertFalse(e.isAnchored(5));
		assertFalse(e.isAnchored(6));
		assertEquals(sce.getEvidenceID(), e.evidence().getEvidenceID());
		for (int i = 0; i < e.length(); i++) {
			assertEquals(2 + i, e.node(i).lastStart());
			assertEquals(2 + i, e.node(i).lastEnd());
			assertEquals(e.isAnchored(i), e.node(i).isReference());
			assertEquals(e.kmer(i), e.node(i).lastKmer());
			assertEquals(e.weight(i), e.node(i).weight());
			assertEquals(e, e.node(i).evidence());
		}
	}
	@Test
	public void bwd_softclip() {
		SoftClipEvidence sce = SCE(BWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1S4M6S"))));
		KmerEvidence e = KmerEvidence.create(4, sce);
		assertEquals(1, e.startPosition());
		assertEquals(1, e.endPosition());
		assertEquals(11-6-(4-1), e.length());
		assertEquals(1, e.weight(0));
		assertEquals(2, e.weight(1));
		assertEquals("ACGT", K(4, e.kmer(0)));
		assertEquals("CGTT", K(4, e.kmer(1)));
		assertFalse(e.isAnchored(0));
		assertTrue(e.isAnchored(1));
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
	private void assertMatches(SingleReadEvidence ex, SingleReadEvidence act) {
		KmerEvidence expected = KmerEvidence.create(4, ex);
		KmerEvidence actual = KmerEvidence.create(4, act);
		assertEquals(expected.length(), actual.length());
		for (int i = 0; i < expected.length(); i++) {
			assertEquals(expected.kmer(i), actual.kmer(i));
			assertEquals(expected.weight(i), actual.weight(i));
			assertEquals(expected.isAnchored(i), actual.isAnchored(i));
		}
	}
	@Test
	public void indel() {
		assertMatches(SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "5M6S")))[0]),
				IE(withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "5M1D6M")))[0]));
		
		assertMatches(SCE(BWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 7, "6S5M")))[0]),
				IE(withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "1M4D5I5M")))[0]).asRemote());
	}
	@Test
	public void split_read() {
		assertMatches(SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "5M6S")))[0]),
				SR(withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 2, "5M6S")))[0], Read(0, 4, "6M")));
		
		assertMatches(SCE(BWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 7, "6S5M")))[0]),
				SR(withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTTATACCG", Read(0, 7, "6S5M")))[0], Read(0, 3, "2S2M2S")));
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
		assertEquals(4, KmerEvidence.create(2, SCE(FWD, Read(0, 1, "10H1S2M3S"))).length());
		assertEquals(2, KmerEvidence.create(2, SCE(BWD, Read(0, 1, "10H1S2M3S"))).length());
		assertEquals(7, KmerEvidence.create(2, SCE(FWD, Read(0, 1, "5M3S"))).length());
		assertEquals(7, KmerEvidence.create(2, SCE(BWD, Read(0, 1, "3S5M"))).length());
		
	}
	@Test
	public void should_return_null_if_kmer_longer_than_read() {
		assertNotNull(KmerEvidence.create(12, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")))));
		assertNull(KmerEvidence.create(13, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")))));
		assertNotNull(KmerEvidence.create(11, NRRP(withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
		assertNull(KmerEvidence.create(12, NRRP(withSequence("ACGTTATACCG", DP(0, 21, "1S9M1S", false, 1, 1, "11M", false)))));
	}
	@Test
	public void short_sc_anchor_should_be_considered_unanchored_evidence() {
		KmerEvidence e = KmerEvidence.create(3, SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "2M10S"))));
		assertTrue(IntStream.range(0, e.length()).allMatch(i -> !e.isAnchored(i)));
	}
	@Test
	public void should_exclude_ambiguous_kmers() {
		KmerEvidence e = KmerEvidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("ACGTNATACCG", Read(0, 2, "5M6S")))));
		assertNotNull(e.node(0));
		assertNull(e.node(1));
		assertNull(e.node(2));
		assertNull(e.node(3));
		assertNull(e.node(4));
		assertNotNull(e.node(5));
		e = KmerEvidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4}, withSequence("ACGTN", Read(0, 2, "4M1S")))));
		assertNotNull(e.node(0));
		assertNull(e.node(1));
		e = KmerEvidence.create(4, SCE(FWD, withQual(new byte[] { 0,1,2,3,4,5,6,7,8,9,10}, withSequence("NCGTTATACCN", Read(0, 2, "5M6S")))));
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
	@Test
	public void isAnchor_should_be_false_if_any_kmer_base_falls_outside_reference_contigs_bounds() {
		KmerEvidence kestart = KmerEvidence.create(2, SCE(FWD, Read(0, -1, "100M5S"))); 
		Assert.assertFalse(kestart.isAnchored(0));
		Assert.assertFalse(kestart.isAnchored(1));
		Assert.assertTrue(kestart.isAnchored(2));
		Assert.assertTrue(kestart.isAnchored(3));
		//    end
		//     |
		// 67890
		//  SMMMMMM
		KmerEvidence keend = KmerEvidence.create(2, SCE(BWD, Read(0, 9998, "1S6M")));
		Assert.assertFalse(keend.isAnchored(0));
		Assert.assertTrue(keend.isAnchored(1));
		Assert.assertTrue(keend.isAnchored(2));
		Assert.assertFalse(keend.isAnchored(3));
		Assert.assertFalse(keend.isAnchored(4));
	}
	@Test
	public void readpair_breakpoint_contig_bounds_truncation_should_not_affect_kmer_placement() {
		//      |<-- start of contig
		//543210
		//      1234567890
		//      MMMM
		//  ********
		//**********
		MockSAMEvidenceSource ses = SES(8, 10);
		KmerEvidence e = KmerEvidence.create(2, NRRP(ses, withSequence("ACGT", DP(0, 1, "4M", false, 1, 1, "4M", true))));
		assertEquals(-5, e.startPosition());
		assertEquals(-3, e.endPosition());
	}
	@Test(expected=NotImplementedException.class)
	public void does_not_handle_XNX_unanchored_clips() {
		KmerEvidence.create(2, SCE(BWD, Read(0, 1, "10S1X10N1X")));
	}
	@Test
	public void flagSelfIntersectingKmersAsAmbiguous_should_flag_kmers_furtherest_from_forward_anchor() {
		MockSAMEvidenceSource ses = SES(0, 100);
		ses.getContext().getAssemblyParameters().positional.trimSelfIntersectingReads = true;
		KmerEvidence e = KmerEvidence.create(4, NRRP(ses, withSequence("AATTAATTAATT", DP(0, 1, "12M", true, 1, 1, "12M", true))));
		// AATT 0 
		//  ATTA 1
		//   TTAA 2
		//    TAAT 3 - loops back to AATT - trim from here onward
		for (int i = 0; i < 3; i++) {
			assertNotNull(e.node(i));
		}
		for (int i = 3; i < e.length(); i++) {
			assertNull(e.node(i));
		}
	}
	@Test
	public void flagSelfIntersectingKmersAsAmbiguous_should_flag_kmers_furtherest_from_backward_anchor() {
		MockSAMEvidenceSource ses = SES(0, 100);
		ses.getContext().getAssemblyParameters().positional.trimSelfIntersectingReads = true;
		KmerEvidence e = KmerEvidence.create(4, NRRP(ses, withSequence("AATTAATTAATT", DP(0, 1, "12M", false, 1, 1, "12M", true))));
		for (int i = e.length(); i > e.length() - 4; i--) {
			assertNotNull(e.node(i));
		}
		for (int i = e.length() - 4; i >=0; i--) {
			assertNull(e.node(i));
		}
	}
	@Test
	public void read_pair_createAnchor_should_place_read_at_mapped_position() {
		assertEquals(100, KmerEvidence.createAnchor(1, NRRP(SES(2, 10), OEA(0, 100, "1M", false)), 0, null).startPosition());
		assertEquals(100, KmerEvidence.createAnchor(1, NRRP(SES(2, 10), OEA(0, 100, "1M", false)), 0, null).endPosition());
	}
	@Test
	public void anchor_should_not_be_equal() {
		NonReferenceReadPair rp = NRRP(SES(), OEA(0, 100, "1M2I3M", false));
		KmerEvidence evidence = KmerEvidence.create(1, rp);
		KmerEvidence anchor = KmerEvidence.createAnchor(1, rp, 0, null);
		assertFalse(evidence.equals(anchor));
		assertNotEquals(evidence.hashCode(), anchor.hashCode());
	}
	@Test
	public void read_pair_should_place_relative_to_anchoring_base_closest_to_break() {
		//           2         3
		//  123456789012345678901234567890
		//         ssMddddddddddM  break is after alignment so we define relative to the final anchoring base at pos 31
		//                   mmmm <- inferred frag from anchor
		//                   |    <- with fragment size of one this puts the fragment here
		//                MMMM    <- and our mate is relative to fragment so it has to be at least 25
		SAMEvidenceSource ses = SES(1, 10);
		KmerEvidence ke = KmerEvidence.create(1, NRRP(ses, OEA(0, 20, "2S1M10D1M", true)));
		assertTrue(ke.node(0).lastStart() >= 25);

	}
	@Test
	public void read_pair_must_have_at_least_one_base_after_anchoring_base_closest_to_break() {
		SAMEvidenceSource ses = SES(1, 10);
		//           2
		//  12345678901234567890
		//         ssMddddM       break occurs after 25 so mate must have at least one base after 26
		//              MMMM
		//                >>>> or further this direction
		KmerEvidence ke = KmerEvidence.create(1, NRRP(ses, OEA(0, 20, "2S1M4D1M", true)));
		assertEquals(23, ke.node(0).lastStart());
		//           2
		//  12345678901234567890
		//           M----Mss
		//          MMMM
		//           <<<< or further this direction
		ke = KmerEvidence.create(1, NRRP(ses, OEA(0, 20, "1M4D1M2S", false)));
		assertEquals(19, ke.node(0).lastEnd());
	}

	/**
	 * Regression edge case in which KmerEvidence was using {start, evidenceID} for equality.
	 * Turns out with just the right length anchoring soft clip, you can get both the anchor
	 * and the read to start at the same position.
	 */
	@Test
	public void rp_should_have_anchor_and_record_not_equal_issue349() {
		NonReferenceReadPair dp = NRRP(SES(0, 1000), DP(0, 100, "73S26M1S", true, 0, 51, "1S99M", true));
		KmerEvidence ke = KmerEvidence.create(2, dp);
		KmerEvidence keAnchor = KmerEvidence.createAnchor(1, dp, 0, SMALL_FA);
		Assert.assertNotSame(ke, keAnchor);
		Assert.assertNotEquals(ke, keAnchor);
	}

	@Test
	public void issue314_rp_bounds_calculation_should_handle_unequal_read_lengths() {
		NonReferenceReadPair dp = NRRP(SES(0, 12), DP(0, 100, "2S2M3S", true, 0, 51, "3M", true));
		KmerEvidence ke = KmerEvidence.create(2, dp);
		//          1
		//          0
		// 12345678901234567890
		//        ssMMsss
		//        |----------| max 12bp fragment
		//          ***        left-most possible alignment (at least one base must not overlap the unclipped alignment)
		//                 *** right-most possible alignment
		Assert.assertEquals(100, ke.node(0).firstStart());
		Assert.assertEquals(107, ke.node(0).firstEnd());

		dp = NRRP(SES(0, 12), DP(0, 100, "2S2M3S", false, 0, 51, "3M", true));
		ke = KmerEvidence.create(2, dp);
		//                    1
		//          9         0
		// 123456789012345678901234567890
		//                  ssMMsss
		//             |----------| max 12bp fragment
		//             ***          left-most possible alignment
		//                   ***    right-most possible alignment
		Assert.assertEquals(93, ke.node(0).firstStart());
		Assert.assertEquals(99, ke.node(0).firstEnd());
	}
}
