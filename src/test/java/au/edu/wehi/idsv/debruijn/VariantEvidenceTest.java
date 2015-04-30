package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.RealignedSoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class VariantEvidenceTest extends TestHelper {
	// FSC:
	// 0        1
	// 1234567890123456789
	// SSSSMMMMMMMSSSSSSSS start=4, ref=7, end=8, k=4
	// SSSS      |       | 0  ignored kmer
	//  SSSM     |       | 1  ignored kmer
	//   SSMM    |       | 2  ignored kmer
	//    SMMM   |       | 3  ignored kmer
	//	   MMMM  |       | 4  ref kmer
	//      MMMM |       | 5  ref kmer
	//       MMMM|       | 6  ref kmer
	//        MMMM       | 7  ref kmer     <- anchor kmer
	//         MMMS      | 8  branch kmer
	//          MMSS     | 9  branch kmer
	//           MSSS    | 10 branch kmer
	//           |SSSS   | 11 branch kmer
	//           | SSSS  | 12 branch kmer
	//           |  SSSS | 13 branch kmer
	//           |   SSSS| 14 branch kmer
	//           |    SSSS 15 branch kmer
	//           ^
	//        anchor pos
	public VariantEvidence fsc = new VariantEvidence(4, SCE(FWD, Read(0, 5, "4S7M8S")), getContext().getLinear());
	public VariantEvidence fdp = new VariantEvidence(4, NRRP(DP(0, 1, "2M", true, 1, 500, "4S7M8S", false)), getContext().getLinear());
	// RSC:
	// 0        1
	// 1234567890123456789
	// SSSSMMMMMMMSSSSSSSS start=4, ref=7, end=8, k=4
	// SSSS|             | 0  branch kmer
	//  SSSM             | 1  branch kmer
	//   SSMM            | 2  branch kmer
	//    SMMM           | 3  branch kmer
	//	   MMMM          | 4  ref kmer     <- anchor kmer
	//     |MMMM         | 5  ref kmer
	//     | MMMM        | 6  ref kmer
	//     |  MMMM       | 7  ref kmer
	//     |   MMMS      | 8  ignored kmer
	//     |    MMSS     | 9  ignored kmer
	//     |     MSSS    | 10 ignored kmer
	//     |      SSSS   | 11 ignored kmer
	//     |       SSSS  | 12 ignored kmer
	//     |        SSSS | 13 ignored kmer
	//     |         SSSS| 14 ignored kmer
	//     |          SSSS 15 ignored kmer
	//     ^
	//  anchor pos
	public VariantEvidence rsc = new VariantEvidence(4, SCE(BWD, Read(0, 5, "4S7M8S")), getContext().getLinear());
	public VariantEvidence rdp = new VariantEvidence(4, NRRP(DP(0, 1, "2M", false, 1, 500, "4S7M8S", true)), getContext().getLinear());
	@Test
	public void getReadKmers_should_read_kmers_in_ascending_genomic_coordinates() {
		VariantEvidence e = new VariantEvidence(2, SCE(FWD, withSequence("ACGTTTTT", Read(0, 5, "4M4S"))), getContext().getLinear());
		List<ReadKmer> kmers = Lists.newArrayList(e.getReadKmers());
		assertEquals(7, kmers.size());
		assertEquals("AC", S(KmerEncodingHelper.encodedToPicardBases(2, kmers.get(0).kmer)));
		assertEquals("CG", S(KmerEncodingHelper.encodedToPicardBases(2, kmers.get(1).kmer)));
	}
	@Test
	public void shouldReverseComp_for_concordant_position() {
		assertEquals("GA", S(KmerEncodingHelper.encodedToPicardBases(2, Lists.newArrayList(new VariantEvidence(2, NRRP(withSequence("TTTC", DP(0, 1, "4M", true, 1, 5, "4M", true))), getContext().getLinear()).getReadKmers()).get(0).kmer)));
		assertEquals("TT", S(KmerEncodingHelper.encodedToPicardBases(2, Lists.newArrayList(new VariantEvidence(2, NRRP(withSequence("TTTC", DP(0, 1, "4M", true, 1, 5, "4M", false))), getContext().getLinear()).getReadKmers()).get(0).kmer)));
		assertEquals("TT", S(KmerEncodingHelper.encodedToPicardBases(2, Lists.newArrayList(new VariantEvidence(2, NRRP(withSequence("TTTC", DP(0, 1, "4M", false, 1, 5, "4M", true))), getContext().getLinear()).getReadKmers()).get(0).kmer)));
		assertEquals("GA", S(KmerEncodingHelper.encodedToPicardBases(2, Lists.newArrayList(new VariantEvidence(2, NRRP(withSequence("TTTC", DP(0, 1, "4M", false, 1, 5, "4M", false))), getContext().getLinear()).getReadKmers()).get(0).kmer)));
	}
	@Test
	public void isReferenceKmer_should_require_all_bases_matching_reference() {
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 4 && i <= 7, fsc.isReferenceKmer(i));
			assertEquals(i >= 4 && i <= 7, rsc.isReferenceKmer(i));
			assertEquals(false, fdp.isReferenceKmer(i));
			assertEquals(false, rdp.isReferenceKmer(i));
		}
	}
	@Test
	public void isSkippedKmer_should_skip_nonvariant_SC_for_SC_evidence() {
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 0 && i <= 3, fsc.isSkippedKmer(i));
			assertEquals(i >= 8 && i <= 15, rsc.isSkippedKmer(i));
		}
	}
	@Test
	public void isSkippedKmer_should_not_skip_read_SC_pair_bases() {
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(false, fdp.isSkippedKmer(i));
			assertEquals(false, rdp.isSkippedKmer(i));
		}
	}
	@Test
	public void getExpectedLinearPosition_should_return_linear_genomic_position_as_if_concordantly_mapped() {
		// 12345
		// **-MM
		VariantEvidence v = new VariantEvidence(4, NRRP(SES(5, 5), DP(0, 4, "2M", false, 1, 500, "2M", true)), getContext().getLinear());
		assertEquals(LCCB + 1, v.getExpectedLinearPosition(0));
		
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(LCCB + 1 + i, fsc.getExpectedLinearPosition(i));
			assertEquals(LCCB + 1 + i, rsc.getExpectedLinearPosition(i));
			assertEquals(LCCB + 1 + 300 - 19 + i, fdp.getExpectedLinearPosition(i));
			assertEquals(LCCB + 1 - 300 + i + 1 + 1, rdp.getExpectedLinearPosition(i));
		}
	}
	@Test
	public void sc_not_dp_should_be_directly_anchored_to_reference() {
		assertFalse(fdp.isDirectlyAnchoredToReference());
		assertFalse(rdp.isDirectlyAnchoredToReference());
		assertTrue(fsc.isDirectlyAnchoredToReference());
		assertTrue(rsc.isDirectlyAnchoredToReference());
	}
	@Test
	public void getReferenceKmerCount_should_count_reference_supporting_kmers() {
		assertEquals(4, fsc.getReferenceKmerCount());
		assertEquals(4, rsc.getReferenceKmerCount());
		assertEquals(0, fdp.getReferenceKmerCount());
		assertEquals(0, rdp.getReferenceKmerCount());
	}
	@Test
	public void getExpectedLinearPosition_should_be_for_local_breakend() {
		for (VariantEvidence v : new VariantEvidence[] { 
				new VariantEvidence(4, NRRP(SES(5, 5), DP(0, 1, "2M", true, 0, 1000, "2M", false)), getContext().getLinear()),
				new VariantEvidence(4, SCE(FWD, Read(0, 5, "5M5S")), getContext().getLinear()),
				new VariantEvidence(4, SCE(FWD, Read(0, 5, "5M5S"), Read(0, 100, "5M")), getContext().getLinear()),
				new VariantEvidence(4, ((RealignedSoftClipEvidence)SCE(BWD, Read(0, 100, "5S5M"), Read(0, 5, "5M"))).asRemote(), getContext().getLinear())
				}) {
			for (int i = 0; i < Lists.newArrayList(v.getReadKmers()).size(); i++) {
				assertTrue(v.getExpectedLinearPosition(i) < LCCB + 350);
				assertTrue(v.getExpectedLinearPosition(i) > LCCB);
			}
		}
	}
	@Test
	public void getExpectedLinearPosition_should_assume_reference_mapping() {
		for (VariantEvidence v : new VariantEvidence[] {
				new VariantEvidence(4, SCE(FWD, Read(0, 1, "10M5S")), getContext().getLinear()),
				new VariantEvidence(4, SCE(BWD, Read(0, 6, "5S10M")), getContext().getLinear()),
				}) {
			for (int i = 0; i < Lists.newArrayList(v.getReadKmers()).size(); i++) {
				assertEquals(LCCB + i + 1, v.getExpectedLinearPosition(i));
			}
		}
	}
}