package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class SAMRecordUtilTest extends TestHelper {
	@Test
	public void isAlignmentSoftClipped() {
		assertTrue(SAMRecordUtil.isAlignmentSoftClipped(Read(0, 1, "1S1M")));
		assertTrue(SAMRecordUtil.isAlignmentSoftClipped(Read(0, 1, "1S1M1S")));
		assertTrue(SAMRecordUtil.isAlignmentSoftClipped(Read(0, 1, "1H1S1M1S1H")));
		assertFalse(SAMRecordUtil.isAlignmentSoftClipped(Read(0, 1, "1H1M")));
		assertFalse(SAMRecordUtil.isAlignmentSoftClipped(Read(0, 1, "1H1M1H")));
		assertFalse(SAMRecordUtil.isAlignmentSoftClipped(Read(0, 1, "1M1D1I1M")));
	}
	@Test
	public void getStartSoftClipLength() {
		assertEquals(0, SAMRecordUtil.getStartSoftClipLength(Read(0, 1, "1M1S")));
		assertEquals(1, SAMRecordUtil.getStartSoftClipLength(Read(0, 1, "1S1M1S")));
		assertEquals(2, SAMRecordUtil.getStartSoftClipLength(Read(0, 1, "2S1M1S")));
		assertEquals(2, SAMRecordUtil.getStartSoftClipLength(Read(0, 1, "1H2S1M1S")));
	}
	@Test
	public void getEndSoftClipLength() {
		assertEquals(0, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1S1M")));
		assertEquals(1, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1S1M1S")));
		assertEquals(2, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1S1M2S")));
		assertEquals(2, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "2M2S1H")));
	}
	@Test
	public void ensureNmTag_should_not_require_reference_if_tag_set() {
		SAMRecord r = new SAMRecord(null);
		r.setAttribute("NM", 1);
		SAMRecordUtil.ensureNmTag((ReferenceSequenceFileWalker)null, r);
	}

	@Test
	public void ensureNmTag_should_set_NM() {
		// AAAA
		// ACGT
		SAMRecord r = Read(1, 1, "4M");
		r.clearAttributes();
		SAMRecordUtil.ensureNmTag(new ReferenceSequenceFileWalker(SMALL_FA), r);
		assertEquals(3, (int)r.getIntegerAttribute("NM"));
	}
	@Test
	public void getStartSoftClipBases() {
		assertEquals("", S(SAMRecordUtil.getStartSoftClipBases(withSequence("ATGC", Read(0, 1, "1H1M3S1H"))[0])));
		assertEquals("A", S(SAMRecordUtil.getStartSoftClipBases(withSequence("ATGC", Read(0, 1, "1S1M2S"))[0])));
		assertEquals("AT", S(SAMRecordUtil.getStartSoftClipBases(withSequence("ATGC", Read(0, 1, "2S1M1S2H"))[0])));
		assertEquals("ATG", S(SAMRecordUtil.getStartSoftClipBases(withSequence("ATGC", Read(0, 1, "3S1M10H"))[0])));
	}
	@Test
	public void getEndSoftClipBases() {
		assertEquals("TGC", S(SAMRecordUtil.getEndSoftClipBases(withSequence("ATGC", Read(0, 1, "1H1M3S1H"))[0])));
		assertEquals("GC", S(SAMRecordUtil.getEndSoftClipBases(withSequence("ATGC", Read(0, 1, "1S1M2S"))[0])));
		assertEquals("C", S(SAMRecordUtil.getEndSoftClipBases(withSequence("ATGC", Read(0, 1, "10H2S1M1S"))[0])));
		assertEquals("", S(SAMRecordUtil.getEndSoftClipBases(withSequence("ATGC", Read(0, 1, "3S1M"))[0])));
	}
	@Test
	public void getStartSoftClipBaseQualities() {
		assertArrayEquals(new byte[] {}, SAMRecordUtil.getStartSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "1H1M3S1H"))[0]));
		assertArrayEquals(new byte[] {1}, SAMRecordUtil.getStartSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "1S1M2S"))[0]));
		assertArrayEquals(new byte[] {1,2}, SAMRecordUtil.getStartSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "10H2S1M1S"))[0]));
		assertArrayEquals(new byte[] {1,2,3}, SAMRecordUtil.getStartSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "3S1M"))[0]));
	}
	@Test
	public void getEndSoftClipBaseQualities() {
		assertArrayEquals(new byte[] {2,3,4}, SAMRecordUtil.getEndSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "1H1M3S1H"))[0]));
		assertArrayEquals(new byte[] {3,4}, SAMRecordUtil.getEndSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "1S1M2S"))[0]));
		assertArrayEquals(new byte[] {4}, SAMRecordUtil.getEndSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "10H2S1M1S"))[0]));
		assertArrayEquals(new byte[] {}, SAMRecordUtil.getEndSoftClipBaseQualities(withQual(new byte[] {1,2,3,4}, Read(0, 1, "3S1M"))[0]));
	}
	@Test
	public void getSoftClipBaseQualities_should_handle_both_nulls() {
		assertNull(SAMRecordUtil.getEndSoftClipBaseQualities(withQual(SAMRecord.NULL_QUALS, Read(0, 1, "1H1M3S1H"))[0]));
		assertNull(SAMRecordUtil.getEndSoftClipBaseQualities(withQual(null, Read(0, 1, "1H1M3S1H"))[0]));
		
		assertNull(SAMRecordUtil.getStartSoftClipBaseQualities(withQual(SAMRecord.NULL_QUALS, Read(0, 1, "1H1M3S1H"))[0]));
		assertNull(SAMRecordUtil.getStartSoftClipBaseQualities(withQual(null, Read(0, 1, "1H1M3S1H"))[0]));
	}
	@Test
	public void getSoftClipBases_should_handle_both_nulls() {
		assertNull(SAMRecordUtil.getEndSoftClipBases(withSequence((byte[])null, Read(0, 1, "1S1M2S"))[0]));
		assertNull(SAMRecordUtil.getEndSoftClipBases(withSequence(SAMRecord.NULL_SEQUENCE, Read(0, 1, "10H2S1M1S"))[0]));
		
		assertNull(SAMRecordUtil.getStartSoftClipBases(withSequence((byte[])null, Read(0, 1, "1S1M2S"))[0]));
		assertNull(SAMRecordUtil.getStartSoftClipBases(withSequence(SAMRecord.NULL_SEQUENCE, Read(0, 1, "10H2S1M1S"))[0]));
	}
	@Test
	public void getTotalReferenceBaseQual_should_sum_mapped_bases() {
		byte[] qual = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		assertEquals(3+4+5+6+7, SAMRecordUtil.getTotalReferenceBaseQual(withQual(qual, Read(0, 1, "3S5M2S"))[0]));
	}
	@Test
	public void getMaxReferenceBaseQual_should_max_mapped_bases() {
		byte[] qual = { 0, 1, 2, 3, 4, 5, 7, 6, 8, 9 };
		assertEquals(7, SAMRecordUtil.getMaxReferenceBaseQual(withQual(qual, Read(0, 1, "3S5M2S"))[0]));
	}
	@Test
	public void overlap() {
		assertTrue(SAMRecordUtil.overlap(DP(1, 1, "5M", true, 1, 5, "5M", false)[0], DP(1, 1, "5M", true, 1, 5, "5M", false)[1]));
		assertTrue(SAMRecordUtil.overlap(DP(1, 1, "5M", true, 1, 5, "5M", true)[0], DP(1, 1, "5M", true, 1, 5, "5M", true)[1]));
		assertFalse(SAMRecordUtil.overlap(DP(1, 1, "5M", true, 0, 5, "5M", false)[0], DP(1, 1, "5M", true, 0, 5, "5M", false)[1]));
		assertFalse(SAMRecordUtil.overlap(DP(1, 1, "5M", true, 1, 6, "5M", true)[0], DP(1, 1, "5M", true, 1, 6, "5M", true)[1]));
		
		assertFalse(SAMRecordUtil.overlap(OEA(0, 1, "5M", true)[0], OEA(0, 1, "5M", true)[1]));
	}
	@Test
	public void estimateFragmentSize() {
		assertEquals(1, SAMRecordUtil.estimateFragmentSize(RP(0, 1, 1, 1)[0], PairOrientation.FR));
		assertEquals(1, SAMRecordUtil.estimateFragmentSize(RP(0, 1, 1, 1)[1], PairOrientation.FR));
		assertEquals(2, SAMRecordUtil.estimateFragmentSize(RP(0, 1, 2, 1)[0], PairOrientation.FR));
		assertEquals(2, SAMRecordUtil.estimateFragmentSize(RP(0, 1, 2, 1)[1], PairOrientation.FR));
		assertEquals(3, SAMRecordUtil.estimateFragmentSize(RP(0, 1, 2, 2)[0], PairOrientation.FR));
		assertEquals(3, SAMRecordUtil.estimateFragmentSize(RP(0, 1, 2, 2)[1], PairOrientation.FR));
	}
	@Test
	public void estimateFragmentSize_should_assume_FR_orientation() {
		SAMRecord[] dp = DP(0, 1, "1M", true, 0, 5, "1M", true);
		assertEquals(0, SAMRecordUtil.estimateFragmentSize(dp[0], PairOrientation.FR));
		assertEquals(0, SAMRecordUtil.estimateFragmentSize(dp[1], PairOrientation.FR));
		dp = DP(0, 1, "1M", false, 0, 5, "1M", false);
		assertEquals(0, SAMRecordUtil.estimateFragmentSize(dp[0], PairOrientation.FR));
		assertEquals(0, SAMRecordUtil.estimateFragmentSize(dp[1], PairOrientation.FR));
	}
	@Test
	public void calculateFragmentSize_should_assume_FR_orientation() {
		SAMRecord[] dp = DP(0, 1, "1M", true, 0, 5, "1M", true);
		assertEquals(0, SAMRecordUtil.calculateFragmentSize(dp[0], dp[1], PairOrientation.FR));
		assertEquals(0, SAMRecordUtil.calculateFragmentSize(dp[1], dp[0], PairOrientation.FR));
		dp = DP(0, 1, "1M", false, 0, 5, "1M", false);
		assertEquals(0, SAMRecordUtil.calculateFragmentSize(dp[0], dp[1], PairOrientation.FR));
		assertEquals(0, SAMRecordUtil.calculateFragmentSize(dp[1], dp[0], PairOrientation.FR));
	}
	@Test
	public void fragmentSize_estimation_should_match_calculation_for_mapped_reads() {
		for (int i = 1; i < 16; i++) {
			SAMRecord[] rp = RP(0, 1, i, 1);
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[0], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[0], rp[1], PairOrientation.FR));
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[1], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[0], rp[1], PairOrientation.FR));
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[0], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[1], rp[0], PairOrientation.FR));
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[1], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[1], rp[0], PairOrientation.FR));
		}
		for (int i = 2; i < 16; i++) {
			SAMRecord[] rp = RP(0, 1, i, 2);
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[0], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[0], rp[1], PairOrientation.FR));
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[1], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[0], rp[1], PairOrientation.FR));
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[0], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[1], rp[0], PairOrientation.FR));
			assertEquals(SAMRecordUtil.estimateFragmentSize(rp[1], PairOrientation.FR), SAMRecordUtil.calculateFragmentSize(rp[1], rp[0], PairOrientation.FR));
		}
	}
	@Test
	public void estimatedReadsOverlap_should_call_conservatively_and_assume_18_mapped_bases_in_mate() {
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 80, 20)[0], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 80, 20)[1], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 81, 20)[0], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 81, 20)[1], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 99, 20)[0], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 99, 20)[1], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 100, 20)[0], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 100, 20)[1], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 110, 20)[0], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 110, 20)[1], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 117, 20)[0], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 117, 20)[1], PairOrientation.FR, 18));
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 118, 20)[0], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 118, 20)[1], PairOrientation.FR, 18)); // conservate call since we don't know mate CIGAR
		assertTrue(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 119, 20)[0], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 119, 20)[1], PairOrientation.FR, 18)); // conservate call since we don't know mate CIGAR
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 120, 20)[0], PairOrientation.FR, 18));
		assertFalse(SAMRecordUtil.estimatedReadsOverlap(RP(0, 100, 120, 20)[1], PairOrientation.FR, 18));
	}
	@Test
	public void isDovetailing_should_require_overlap() {
		assertFalse(SAMRecordUtil.isDovetailing(RP(0, 1, 2, 1)[0], RP(0, 1, 2, 1)[1], PairOrientation.FR, 4));
	}
	@Test
	public void trimSoftClips_should_remove_from_both_end() {
		SAMRecord r = Read(0,  1, "3S1M1S");
		r.setReadBases(B("ACGTN"));
		r.setBaseQualities(B("ACGTN"));
		SAMRecordUtil.trimSoftClips(r, 2, 1);
		assertEquals("1S1M", r.getCigarString());
		assertEquals("GT", S(r.getReadBases()));
		assertEquals("GT", S(r.getBaseQualities()));
	}
	@Test
	public void isReferenceAlignment_should_allow_only_match_mismatch_and_hard_clipping() {
		assertTrue(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1M")));
		assertTrue(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1H1X1=")));
		assertFalse(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1S1M")));
		assertFalse(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1M1N1M")));
		assertFalse(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1M1D1M")));
		assertFalse(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1M1I1M")));
	}
	@Test
	public void isReferenceAlignment_should_ignore_zero_length_elements() {
		assertTrue(SAMRecordUtil.isReferenceAlignment(Read(0, 1, "1M0S")));
	}
	@Test
	public void alignedEntropy_should_return_entropy_excluding_soft_clipped_bases() {
		assertEquals(0, SAMRecordUtil.alignedEntropy(withSequence("ATG", Read(0, 1, "1S1M1S"))[0]), 0);
		assertEquals(0, SAMRecordUtil.alignedEntropy(withSequence("ATTG", Read(0, 1, "1S2M1S"))[0]), 0);
		assertEquals(2, SAMRecordUtil.alignedEntropy(withSequence("AACGTG", Read(0, 1, "1S4M1S"))[0]), 0);
		assertEquals(2, SAMRecordUtil.alignedEntropy(withSequence("AACGTGG", Read(0, 1, "1S4M2S"))[0]), 0);
		assertEquals(2, SAMRecordUtil.alignedEntropy(withSequence("TAACGTGG", Read(0, 1, "2S4M2S"))[0]), 0);
		assertEquals(1, SAMRecordUtil.alignedEntropy(withSequence("TAGTGG", Read(0, 1, "2S2M2S"))[0]), 0);
		assertEquals(2, SAMRecordUtil.alignedEntropy(withSequence("ACGT", Read(0, 1, "4M"))[0]), 0);
	}
	@Test
	public void entropy_should_ignore_cigar_and_alignment() {
		assertEquals(2, SAMRecordUtil.entropy(withSequence("ACGT", Read(0, 1, "3S1M"))[0]), 0);
	}
	@Test
	public void splitAfter_should_make_two_full_length_reads() {
		Pair<SAMRecord, SAMRecord> p = SAMRecordUtil.splitAfter(Read(0, 1, "10M"), 1);
		assertEquals(1, p.getLeft().getAlignmentStart());
		assertEquals("2M8S", p.getLeft().getCigarString());
		assertEquals(10, p.getLeft().getReadLength());
		
		assertEquals(3, p.getRight().getAlignmentStart());
		assertEquals("2S8M", p.getRight().getCigarString());
		assertEquals(10, p.getRight().getReadLength());
	}
	@Test
	public void splitAfter_should_clean_break() {
		Pair<SAMRecord, SAMRecord> p = SAMRecordUtil.splitAfter(Read(0, 1, "4M4I4M"), 6);
		assertEquals(1, p.getLeft().getAlignmentStart());
		assertEquals("4M8S", p.getLeft().getCigarString());
		assertEquals(12, p.getLeft().getReadLength());
		
		assertEquals(5, p.getRight().getAlignmentStart());
		assertEquals("8S4M", p.getRight().getCigarString());
		assertEquals(12, p.getRight().getReadLength());
	}
	@Test
	public void splitAfter_should_not_alter_base_alignments() {
		Pair<SAMRecord, SAMRecord> p = SAMRecordUtil.splitAfter(Read(0, 1, "4M4D4M"), 3);
		assertEquals(1, p.getLeft().getAlignmentStart());
		assertEquals("4M4S", p.getLeft().getCigarString());
		
		assertEquals(9, p.getRight().getAlignmentStart());
		assertEquals("4S4M", p.getRight().getCigarString());
		
		p = SAMRecordUtil.splitAfter(Read(0, 1, "4M4I4M"), 6);
		assertEquals(1, p.getLeft().getAlignmentStart());
		assertEquals(5, p.getRight().getAlignmentStart());
	}
	@Test
	public void realign_should_use_window_around_read_alignment() {
		SAMRecord read = Read(2, 101, "100M");
		read.setReadBases(B(S(RANDOM).substring(100, 199) + "N"));
		SAMRecord realigned = SAMRecordUtil.realign(SMALL_FA, read, 0, true);
		assertTrue(realigned != read);
		assertEquals("99M1S", realigned.getCigarString());
		assertEquals(101, realigned.getAlignmentStart());
	}
}
