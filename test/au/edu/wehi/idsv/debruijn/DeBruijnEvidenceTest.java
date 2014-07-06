package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class DeBruijnEvidenceTest extends TestHelper {
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
	public DeBruijnEvidence fsc = DeBruijnEvidence.createSoftClipEvidence(FWD, 4, SCE(FWD, Read(0, 5, "4S7M8S")));
	public DeBruijnEvidence fdp = DeBruijnEvidence.createRemoteReadEvidence(FWD, 4, NRRP(DP(0, 1, "2M", true, 1, 5, "4S7M8S", false)));
	// RSC:
	// 0        1
	// 1234567890123456789
	// SSSSMMMMMMMSSSSSSSS start=4, ref=7, end=8, k=4
	// SSSS|             | 15 branch kmer
	//  SSSM             | 14 branch kmer
	//   SSMM            | 13 branch kmer
	//    SMMM           | 12 branch kmer
	//	   MMMM          | 11 ref kmer     <- anchor kmer
	//     |MMMM         | 10 ref kmer
	//     | MMMM        | 9  ref kmer
	//     |  MMMM       | 8  ref kmer
	//     |   MMMS      | 7  ignored kmer
	//     |    MMSS     | 6  ignored kmer
	//     |     MSSS    | 5  ignored kmer
	//     |      SSSS   | 4  ignored kmer
	//     |       SSSS  | 3  ignored kmer
	//     |        SSSS | 2  ignored kmer
	//     |         SSSS| 1  ignored kmer
	//     |          SSSS 0  ignored kmer
	//     ^
	//  anchor pos
	public DeBruijnEvidence rsc = DeBruijnEvidence.createSoftClipEvidence(BWD, 4, SCE(BWD, Read(0, 5, "4S7M8S")));
	public DeBruijnEvidence rdp = DeBruijnEvidence.createRemoteReadEvidence(BWD, 4, NRRP(DP(0, 1, "2M", false, 1, 5, "4S7M8S", true)));
	@Test
	public void shouldReverse_backward_direction_so_all_graph_breakends_are_forward() {
		assertEquals(false, fsc.isReversed());
		assertEquals(true, rsc.isReversed());
		assertEquals(false, fsc.isComplemented());
		assertEquals(false, rsc.isComplemented());
	}
	@Test
	public void should_conditionally_reverse__and_complement_remote_reads_so_they_assemble_with_soft_clips() {
		// In expected orientation => no change
		// MMMM>    <----
		assertEquals(false, DeBruijnEvidence.createRemoteReadEvidence(FWD, 4, NRRP(DP(0, 1, "2M", true, 1, 5, "4S7M8S", false))).isReversed());
		assertEquals(false, DeBruijnEvidence.createRemoteReadEvidence(FWD, 4, NRRP(DP(0, 1, "2M", true, 1, 5, "4S7M8S", false))).isComplemented());
		
		// Aligned to wrong strand -> reverse comp
		// MMMM>    ---->
		assertEquals(true, DeBruijnEvidence.createRemoteReadEvidence(FWD, 4, NRRP(DP(0, 1, "2M", true, 1, 5, "4S7M8S", true))).isReversed());
		assertEquals(true, DeBruijnEvidence.createRemoteReadEvidence(FWD, 4, NRRP(DP(0, 1, "2M", true, 1, 5, "4S7M8S", true))).isComplemented());
		
		// In expected orientation, but we're going backward = only reverse
		// ---->    <MMMM
		assertEquals(true, DeBruijnEvidence.createRemoteReadEvidence(BWD, 4, NRRP(DP(0, 1, "2M", false, 1, 5, "4S7M8S", true))).isReversed());
		assertEquals(false, DeBruijnEvidence.createRemoteReadEvidence(BWD, 4, NRRP(DP(0, 1, "2M", false, 1, 5, "4S7M8S", true))).isComplemented());
		
		// Wrong orientation, and going backward = reverse twice & comp
		// <----    <MMMM
		assertEquals(false, DeBruijnEvidence.createRemoteReadEvidence(BWD, 4, NRRP(DP(0, 1, "2M", false, 1, 5, "4S7M8S", false))).isReversed());
		assertEquals(true, DeBruijnEvidence.createRemoteReadEvidence(BWD, 4, NRRP(DP(0, 1, "2M", false, 1, 5, "4S7M8S", false))).isComplemented());
	}
	@Test
	public void basesSupportingReference_should_count_all_ref_bases() {
		assertEquals(0, fsc.basesSupportingReference(0));
		assertEquals(1, fsc.basesSupportingReference(1));
		assertEquals(2, fsc.basesSupportingReference(2));
		assertEquals(3, fsc.basesSupportingReference(3));
		assertEquals(4, fsc.basesSupportingReference(4));
		assertEquals(4, fsc.basesSupportingReference(5));
		assertEquals(4, fsc.basesSupportingReference(6));
		assertEquals(4, fsc.basesSupportingReference(7));
		assertEquals(3, fsc.basesSupportingReference(8));
		assertEquals(2, fsc.basesSupportingReference(9));
		assertEquals(1, fsc.basesSupportingReference(10));
		assertEquals(0, fsc.basesSupportingReference(11));
		assertEquals(0, fsc.basesSupportingReference(12));
		assertEquals(0, fsc.basesSupportingReference(13));
		assertEquals(0, fsc.basesSupportingReference(14));
		assertEquals(0, fsc.basesSupportingReference(15));
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(0, fdp.basesSupportingReference(i));
		}
		assertEquals(0, rsc.basesSupportingReference(0));
		assertEquals(0, rsc.basesSupportingReference(1));
		assertEquals(0, rsc.basesSupportingReference(2));
		assertEquals(0, rsc.basesSupportingReference(3));
		assertEquals(0, rsc.basesSupportingReference(4));
		assertEquals(1, rsc.basesSupportingReference(5));
		assertEquals(2, rsc.basesSupportingReference(6));
		assertEquals(3, rsc.basesSupportingReference(7));
		assertEquals(4, rsc.basesSupportingReference(8));
		assertEquals(4, rsc.basesSupportingReference(9));
		assertEquals(4, rsc.basesSupportingReference(10));
		assertEquals(4, rsc.basesSupportingReference(11));
		assertEquals(3, rsc.basesSupportingReference(12));
		assertEquals(2, rsc.basesSupportingReference(13));
		assertEquals(1, rsc.basesSupportingReference(14));
		assertEquals(0, rsc.basesSupportingReference(15));
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(0, rdp.basesSupportingReference(i));
		}
	}
	@Test
	public void isReferenceKmer_should_require_all_bases_matching_reference() {
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 4 && i <= 7, fsc.isReferenceKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(false, fdp.isReferenceKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 8 && i <= 11, rsc.isReferenceKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(false, rdp.isReferenceKmer(i));
		}
	}
	@Test
	public void isReferenceKmer_should_allow_zero_skipped_bases() {
		assertTrue(DeBruijnEvidence.createSoftClipEvidence(BWD, 3, SCE(BWD, withSequence("TTATG", Read(0, 10, "2S3M")))).isReferenceKmer(0));
	}
	@Test
	public void isSkippedKmer_should_skip_if_any_SC() {
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 0 && i <= 3, fsc.isSkippedKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(false, fdp.isSkippedKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 0 && i <= 7, rsc.isSkippedKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(false, rdp.isSkippedKmer(i));
		}
	}
	@Test
	public void isVariantKmer_if_any_SC() {
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 8 && i <= 15, fsc.isVariantKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(true, fdp.isVariantKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(i >= 12 && i <= 15, rsc.isVariantKmer(i));
		}
		for (int i= 0; i < 19 - 4; i++) {
			assertEquals(true, rdp.isVariantKmer(i));
		}
	}
	@Test
	public void getReferenceStartingPosition_should_return_genomic_position_of_first_base_of_kmer() {
		assertEquals(1, fsc.getReferenceStartingPosition(0));
		assertEquals(2, fsc.getReferenceStartingPosition(1));
		assertEquals(3, fsc.getReferenceStartingPosition(2));
		assertEquals(4, fsc.getReferenceStartingPosition(3));
		assertEquals(5, fsc.getReferenceStartingPosition(4));
		assertEquals(6, fsc.getReferenceStartingPosition(5));
		assertEquals(7, fsc.getReferenceStartingPosition(6));
		assertEquals(8, fsc.getReferenceStartingPosition(7));
		assertEquals(9, fsc.getReferenceStartingPosition(8));
		assertEquals(10, fsc.getReferenceStartingPosition(9));
		assertEquals(11, fsc.getReferenceStartingPosition(10));
		
		assertEquals(19, rsc.getReferenceStartingPosition(0));
		assertEquals(18, rsc.getReferenceStartingPosition(1));
		assertEquals(17, rsc.getReferenceStartingPosition(2));
		assertEquals(16, rsc.getReferenceStartingPosition(3));
		assertEquals(15, rsc.getReferenceStartingPosition(4));
		assertEquals(14, rsc.getReferenceStartingPosition(5));
		assertEquals(13, rsc.getReferenceStartingPosition(6));
		assertEquals(12, rsc.getReferenceStartingPosition(7));
		assertEquals(11, rsc.getReferenceStartingPosition(8));
		assertEquals(10, rsc.getReferenceStartingPosition(9));
		assertEquals(9, rsc.getReferenceStartingPosition(10));
		assertEquals(8, rsc.getReferenceStartingPosition(11));
		assertEquals(7, rsc.getReferenceStartingPosition(12));
		assertEquals(6, rsc.getReferenceStartingPosition(13));
		assertEquals(5, rsc.getReferenceStartingPosition(14));
	}
	@Test
	public void first_last_methods_should_return_zero_based_kmer_index() {
		assertEquals(0, fsc.firstSkippedKmerOffset());
		assertEquals(3, fsc.lastSkippedKmerOffset());
		assertEquals(4, fsc.firstReferenceKmerOffset());
		assertEquals(7, fsc.lastReferenceKmerOffset());
		assertEquals(8, fsc.firstVariantKmerOffset());
		assertEquals(15, fsc.lastVariantKmerOffset());
		assertEquals(15, fsc.lastKmerOffset());
		assertEquals(16, fsc.kmerCount());
		
		assertEquals(-1, fdp.firstSkippedKmerOffset());
		assertEquals(-1, fdp.lastSkippedKmerOffset());
		assertEquals(-1, fdp.firstReferenceKmerOffset());
		assertEquals(-1, fdp.lastReferenceKmerOffset());
		assertEquals(0, fdp.firstVariantKmerOffset());
		assertEquals(15, fdp.lastVariantKmerOffset());
		assertEquals(15, fdp.lastKmerOffset());
		assertEquals(16, fdp.kmerCount());
		
		assertEquals(0, rsc.firstSkippedKmerOffset());
		assertEquals(7, rsc.lastSkippedKmerOffset());
		assertEquals(8, rsc.firstReferenceKmerOffset());
		assertEquals(11, rsc.lastReferenceKmerOffset());
		assertEquals(12, rsc.firstVariantKmerOffset());
		assertEquals(15, rsc.lastVariantKmerOffset());
		assertEquals(15, rsc.lastKmerOffset());
		assertEquals(16, rsc.kmerCount());
		
		assertEquals(-1, rdp.firstSkippedKmerOffset());
		assertEquals(-1, rdp.lastSkippedKmerOffset());
		assertEquals(-1, rdp.firstReferenceKmerOffset());
		assertEquals(-1, rdp.lastReferenceKmerOffset());
		assertEquals(0, rdp.firstVariantKmerOffset());
		assertEquals(15, rdp.lastVariantKmerOffset());
		assertEquals(15, rdp.lastKmerOffset());
		assertEquals(16, rdp.kmerCount());
	}
	@Test
	public void getInferredReferencePosition_should_assume_direct_mapping_to_reference() {
		assertEquals(8, fsc.getInferredReferencePosition(4));
		assertEquals(9, fsc.getInferredReferencePosition(5));
		assertEquals(10, fsc.getInferredReferencePosition(6));
		assertEquals(11, fsc.getInferredReferencePosition(7));
		
		assertEquals(8, rsc.getInferredReferencePosition(8));
		assertEquals(7, rsc.getInferredReferencePosition(9));
		assertEquals(6, rsc.getInferredReferencePosition(10));
		assertEquals(5, rsc.getInferredReferencePosition(11));
	}
}