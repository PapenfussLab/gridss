package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;


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
	public void getStartClipLength() {
		assertEquals(0, SAMRecordUtil.getStartClipLength(Read(0, 1, "1M1S")));
		assertEquals(1, SAMRecordUtil.getStartClipLength(Read(0, 1, "1S1M1S")));
		assertEquals(2, SAMRecordUtil.getStartClipLength(Read(0, 1, "2S1M1S")));
		assertEquals(3, SAMRecordUtil.getStartClipLength(Read(0, 1, "1H2S1M1S")));
	}
	@Test
	public void getEndSoftClipLength() {
		assertEquals(0, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1S1M")));
		assertEquals(1, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1S1M1S")));
		assertEquals(2, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1S1M2S")));
		assertEquals(2, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "2M2S1H")));
	}
	@Test
	public void getEndClipLength() {
		assertEquals(0, SAMRecordUtil.getEndClipLength(Read(0, 1, "1S1M")));
		assertEquals(1, SAMRecordUtil.getEndClipLength(Read(0, 1, "1S1M1S")));
		assertEquals(2, SAMRecordUtil.getEndClipLength(Read(0, 1, "1S1M2S")));
		assertEquals(3, SAMRecordUtil.getEndClipLength(Read(0, 1, "2M2S1H")));
	}
	@Test
	public void ensureNmTag_should_not_require_reference_if_tag_set() {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
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
	public void getAlignedPercentIdentity_should_match_only_mapped_bases() {
		assertEquals(1, SAMRecordUtil.getAlignedIdentity(withNM(withSequence("NAAAAA", Read(0, 1, "1S5M")))[0]), 0);
		assertEquals(1, SAMRecordUtil.getAlignedIdentity(withNM(withSequence("NTAAAT", Read(0, 1, "2S3M1S")))[0]), 0);
		assertEquals(1, SAMRecordUtil.getAlignedIdentity(withNM(withSequence("NATTTA", Read(0, 1, "1S1M3I1M")))[0]), 0);
		assertEquals(1, SAMRecordUtil.getAlignedIdentity(withNM(withSequence("NAGTAC", Read(1, 1, "1S1M1D4M")))[0]), 0);
		assertEquals(0.5, SAMRecordUtil.getAlignedIdentity(withNM(withSequence("NAATT", Read(0, 1, "1S4M")))[0]), 0);
		assertEquals(0, SAMRecordUtil.getAlignedIdentity(withNM(withSequence("ACCCCA", Read(0, 1, "1S4M1S")))[0]), 0);
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
	@Test
	public void softenHardClips_should_extend_bases() {
		SAMRecord read =  Read(0, 1, "5M5S");
		read.setReadBases(B("AACCGACGTA"));
		read.setBaseQualityString("1234567890");
		SAMRecord read2 =  Read(0, 1, "4H1S5M");
		read2.setReadBases(B("GACGTA"));
		read2.setBaseQualityString("567890");
		SAMRecordUtil.softenHardClips(ImmutableList.of(read, read2));
		assertEquals("5M5S", read.getCigarString());
		assertEquals("AACCGACGTA", read.getReadString());
		assertEquals("1234567890", read.getBaseQualityString());
		assertEquals("5S5M", read2.getCigarString());
		assertEquals("AACCGACGTA", read2.getReadString());
		assertEquals("1234567890", read2.getBaseQualityString());
	}
	@Test
	public void softenHardClips_should_consider_strand() {
		SAMRecord read =  Read(0, 1, "5M5S");
		read.setReadBases(B("AACCGACGTA"));
		read.setBaseQualityString("1234567890");
		SAMRecord read3 =  Read(0, 1, "5H5M");
		read3.setReadNegativeStrandFlag(true);
		read3.setReadBases(B("ACGTA"));
		read3.setBaseQualityString("67890");
		SAMRecordUtil.softenHardClips(ImmutableList.of(read, read3));
		assertEquals("TACGTCGGTT", read3.getReadString());
		assertEquals("0987654321", read3.getBaseQualityString());
	}
	@Test
	public void getSegmentIndex_should_use_zero_based_FI() {
		SAMRecord read =  Read(0, 1, "5M");
		read.setAttribute("FI", 5);
		assertEquals(5, SAMRecordUtil.getSegmentIndex(read));
	}
	@Test
	public void getSegmentIndex_should_use_first_of_pair_flag() {
		SAMRecord read =  Read(0, 1, "5M");
		read.setReadPairedFlag(true);
		read.setFirstOfPairFlag(true);
		assertEquals(0, SAMRecordUtil.getSegmentIndex(read));
	}
	@Test
	public void getSegmentIndex_should_use_second_of_pair_flag() {
		SAMRecord read =  Read(0, 1, "5M");
		read.setReadPairedFlag(true);
		read.setSecondOfPairFlag(true);
		assertEquals(1, SAMRecordUtil.getSegmentIndex(read));
	}
	@Test
	public void getSegmentIndex_shoulddefault_to_zero() {
		SAMRecord read =  Read(0, 1, "5M");
		assertEquals(0, SAMRecordUtil.getSegmentIndex(read));
	}
	@Test
	public void getFirstAlignedBaseReadOffset_should_consider_strand() {
		// forward
		assertEquals(0, SAMRecordUtil.getFirstAlignedBaseReadOffset(Read(0, 1, "5M")));
		assertEquals(0, SAMRecordUtil.getFirstAlignedBaseReadOffset(Read(0, 1, "5M10S")));
		assertEquals(1, SAMRecordUtil.getFirstAlignedBaseReadOffset(Read(0, 1, "1H5M10S")));
		assertEquals(3, SAMRecordUtil.getFirstAlignedBaseReadOffset(Read(0, 1, "1H2S5M10S")));
		// reverse
		assertEquals(0, SAMRecordUtil.getFirstAlignedBaseReadOffset(onNegative(Read(0, 1, "5M"))[0]));
		assertEquals(1, SAMRecordUtil.getFirstAlignedBaseReadOffset(onNegative(Read(0, 1, "5M1S"))[0]));
		assertEquals(3, SAMRecordUtil.getFirstAlignedBaseReadOffset(onNegative(Read(0, 1, "10S5M1S2H"))[0]));
	}
	@Test
	public void getLastAlignedBaseReadOffset_should_consider_strand() {
		// forward
		assertEquals(4, SAMRecordUtil.getLastAlignedBaseReadOffset(Read(0, 1, "5M")));
		assertEquals(4, SAMRecordUtil.getLastAlignedBaseReadOffset(Read(0, 1, "5M10S")));
		assertEquals(1, SAMRecordUtil.getLastAlignedBaseReadOffset(Read(0, 1, "1H1M2S")));
		// reverse
		assertEquals(4, SAMRecordUtil.getLastAlignedBaseReadOffset(onNegative(Read(0, 1, "5M"))[0]));
		assertEquals(5, SAMRecordUtil.getLastAlignedBaseReadOffset(onNegative(Read(0, 1, "5M1S"))[0]));
		assertEquals(7, SAMRecordUtil.getLastAlignedBaseReadOffset(onNegative(Read(0, 1, "10S5M1S2H"))[0]));
	}
	@Test
	public void calculateTemplateTags_should_write_FI_TC() {
		SAMRecord read0 = Read(0, 1, "5M");
		read0.setReadPairedFlag(true);
		read0.setFirstOfPairFlag(true);
		
		SAMRecord read1 = Read(0, 1, "5M");
		read1.setReadPairedFlag(true);
		read1.setSecondOfPairFlag(true);
		
		SAMRecord read5 = Read(0, 1, "5M");
		read5.setAttribute("FI", 5);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1, read5), ImmutableSet.of(SAMTag.FI, SAMTag.TC), false);
		
		assertEquals(0, (int)read0.getIntegerAttribute("FI"));
		assertEquals(1, (int)read1.getIntegerAttribute("FI"));
		assertEquals(5, (int)read5.getIntegerAttribute("FI"));
		
		assertEquals(6, (int)read0.getIntegerAttribute("TC"));
		assertEquals(6, (int)read1.getIntegerAttribute("TC"));
		assertEquals(6, (int)read5.getIntegerAttribute("TC"));
	}
	@Test
	public void calculateTemplateTags_should_write_R2_Q2_using_next_segment_modulo_arithmetic() {
		SAMRecord read0 = Read(0, 1, "3M");
		read0.setReadPairedFlag(true);
		read0.setFirstOfPairFlag(true);
		read0.setReadBases(B("AAC"));
		read0.setBaseQualityString("123");
		
		SAMRecord read1 = Read(0, 1, "3M");
		read1.setReadPairedFlag(true);
		read1.setSecondOfPairFlag(true);
		read1.setReadNegativeStrandFlag(true);
		read1.setReadBases(B("ACT"));
		read1.setBaseQualityString("456");
		
		SAMRecord read5 = Read(0, 1, "3M");
		read5.setAttribute("FI", 5);
		read5.setReadBases(B("GGT"));
		read5.setBaseQualityString("890");
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1, read5), ImmutableSet.of(SAMTag.R2, SAMTag.Q2), false);
		
		assertEquals("AGT", read0.getStringAttribute("R2"));
		assertEquals(null, read1.getStringAttribute("R2"));
		assertEquals("AAC", read5.getStringAttribute("R2"));
		
		assertEquals("654", read0.getStringAttribute("Q2"));
		assertEquals(null, read1.getStringAttribute("Q2"));
		assertEquals("123", read5.getStringAttribute("Q2"));
	}
	@Test
	public void calculateTemplateTags_should_write_R2_Q2_based_on_sequencing_order() {
		SAMRecord read0 = Read(0, 1, "3M");
		read0.setReadPairedFlag(true);
		read0.setFirstOfPairFlag(true);
		read0.setReadBases(B("AAC"));
		read0.setBaseQualityString("123");
		
		SAMRecord read1 = Read(1, 1, "3M");
		read1.setReadPairedFlag(true);
		read1.setSecondOfPairFlag(true);
		read1.setReadNegativeStrandFlag(true);
		read1.setReadBases(B("ACT"));
		read1.setBaseQualityString("456");
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(SAMTag.R2, SAMTag.Q2), false);
		
		assertEquals("AGT", read0.getStringAttribute("R2"));
		assertEquals("AAC", read1.getStringAttribute("R2"));
		
		assertEquals("654", read0.getStringAttribute("Q2"));
		assertEquals("123", read1.getStringAttribute("Q2"));
	}
	@Test
	public void calculateTemplateTags_SA_should_ignore_missing_NM() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		
		read0.setAttribute("NM", null);
		read1.setAttribute("NM", null);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(SAMTag.SA), false);
		
		assertEquals("polyA,2,+,1H1M4H,11,", read0.getStringAttribute("SA"));
		assertEquals("polyA,1,+,1M5H,10,", read1.getStringAttribute("SA"));
	}
	@Test
	public void calculateTemplateTags_SA_should_not_be_written_for_nonchimeric_reads() {
		SAMRecord read0 = Read(0, 1, "3M");
		read0.setReadPairedFlag(true);
		read0.setFirstOfPairFlag(true);
		read0.setReadBases(B("AAC"));
		read0.setBaseQualityString("123");
		
		SAMRecord read1 = Read(1, 1, "3M");
		read1.setReadPairedFlag(true);
		read1.setSecondOfPairFlag(true);
		read1.setReadNegativeStrandFlag(true);
		read1.setReadBases(B("ACT"));
		read1.setBaseQualityString("456");
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(SAMTag.SA), false);
		
		assertEquals(null, read0.getStringAttribute("SA"));
		assertEquals(null, read1.getStringAttribute("SA"));
	}
	/**
	 * Conventionally, at a supplementary line, the first element points to the primary line
	 */
	@Test
	public void calculateTemplateTags_SA_should_write_supplementary_alignments_last_for_chimeric_alignments() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		// 1BP gap
		SAMRecord read2 = withMapq(12, withSequence("G", Read(0, 4, "3H1M2H")))[0];
		SAMRecord read3 = withMapq(13, withSequence("ACGTAA", Read(0, 5, "4S2M")))[0];
		SAMRecord alt = withMapq(14, withSequence("ACGTAA", Read(0, 1, "6M")))[0];
		
		read0.setSupplementaryAlignmentFlag(true);
		read1.setSupplementaryAlignmentFlag(true);
		read2.setSupplementaryAlignmentFlag(false); // read 2 is the 'primary' chimeric alignment record
		read3.setSupplementaryAlignmentFlag(true);
		
		read0.setNotPrimaryAlignmentFlag(true);
		read1.setNotPrimaryAlignmentFlag(true);
		read2.setNotPrimaryAlignmentFlag(true);
		read3.setNotPrimaryAlignmentFlag(true);
		
		read0.setAttribute("NM", 0);
		read1.setAttribute("NM", 1);
		read2.setAttribute("NM", 2);
		read3.setAttribute("NM", 3);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1, read2, read3, alt), ImmutableSet.of(SAMTag.SA), false);
		
		assertEquals(null, alt.getStringAttribute("SA"));
		// "0,4,+,3H1M2H,12,2,;0,1,+,1M5H,10,0,;0,2,+,1H1M4H,11,1,;0,5,+,4S2M,13,3,"
		//   read2                  read0            read1          read3
		assertEquals("polyA,4,+,3H1M2H,12,2;polyA,2,+,1H1M4H,11,1;polyA,5,+,4S2M,13,3", read0.getStringAttribute("SA"));
		assertEquals("polyA,4,+,3H1M2H,12,2;polyA,1,+,1M5H,10,0;polyA,5,+,4S2M,13,3", read1.getStringAttribute("SA"));
		assertEquals("polyA,1,+,1M5H,10,0;polyA,2,+,1H1M4H,11,1;polyA,5,+,4S2M,13,3", read2.getStringAttribute("SA"));
		assertEquals("polyA,4,+,3H1M2H,12,2;polyA,1,+,1M5H,10,0;polyA,2,+,1H1M4H,11,1", read3.getStringAttribute("SA"));
	}
}
