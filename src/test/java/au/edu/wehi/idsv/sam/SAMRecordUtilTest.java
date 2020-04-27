package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.util.UngroupingIterator;
import gridss.ComputeSamTags;
import gridss.SanityCheckEvidence;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import org.junit.experimental.categories.Category;


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
		SAMRecordUtil.ensureNmTag(null, r);
	}

	@Test
	public void ensureNmTag_should_set_NM() {
		// AAAA
		// ACGT
		assertEquals(3, (int)SAMRecordUtil.ensureNmTag(SMALL_FA, withAttr("NM", null, Read(1, 1, "4M"))[0]).getIntegerAttribute("NM"));
		assertEquals(0, (int)SAMRecordUtil.ensureNmTag(SMALL_FA, withSequence("ACGT", withAttr("NM", null, Read(1, 1, "4M")))[0]).getIntegerAttribute("NM"));
		assertEquals(1, (int)SAMRecordUtil.ensureNmTag(SMALL_FA, withSequence("CCGT", withAttr("NM", null, Read(1, 1, "4M")))[0]).getIntegerAttribute("NM"));
		assertEquals(1, (int)SAMRecordUtil.ensureNmTag(SMALL_FA, withSequence("GCGT", withAttr("NM", null, Read(1, 1, "4M")))[0]).getIntegerAttribute("NM"));
		assertEquals(1, (int)SAMRecordUtil.ensureNmTag(SMALL_FA, withSequence("TCGT", withAttr("NM", null, Read(1, 1, "4M")))[0]).getIntegerAttribute("NM"));
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
	public void getEndSoftClipLength_should_not_die_on_fully_SC_reads() {
		assertEquals(100, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "100S")));
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
	public void estimateFragmentSize_should_consider_hard_clips() {
		SAMRecord[] dp = DP(0, 10, "1H2M", true, 0, 13, "3M", false);
		// 1
		// 012345678901234567890
		// MM MM
		assertEquals(7, SAMRecordUtil.estimateFragmentSize(dp[0], PairOrientation.FR));
	}
	@Test
	public void estimateFragmentSize_should_use_mate_cigar_if_present() {
		SAMRecord[] dp = DP(0, 10, "2M", true, 0, 13, "3M", false);
		// 1
		// 012345678901234567890
		// MM MMM
		assertEquals(6, SAMRecordUtil.estimateFragmentSize(dp[0], PairOrientation.FR));
		assertEquals(6, SAMRecordUtil.estimateFragmentSize(dp[1], PairOrientation.FR));
	}
	@Test
	public void estimateFragmentSize_should_use_cigar_based_estimation() {
		SAMRecord[] dp = DP(0, 10, "2M", true, 0, 13, "2M", false);
		dp[0].setReadBases(SAMRecord.NULL_SEQUENCE);
		dp[1].setReadBases(SAMRecord.NULL_SEQUENCE);
		dp[0].setAttribute("MC", null);
		dp[1].setAttribute("MC", null);
		assertEquals(5, SAMRecordUtil.estimateFragmentSize(dp[0], PairOrientation.FR));
		assertEquals(5, SAMRecordUtil.estimateFragmentSize(dp[1], PairOrientation.FR));
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
	private void checkDovetail(boolean expected, int margin, SAMRecord... reads) {
		assert(reads.length == 2);
		assertEquals(expected, SAMRecordUtil.isDovetailing(reads[0], reads[1], PairOrientation.FR, margin));
		assertEquals(expected, SAMRecordUtil.isDovetailing(reads[0], PairOrientation.FR, margin));
	}
	@Test
	public void isDovetailing_should_require_overlap() {
		checkDovetail(true, 4, RP(0, 1, 1, 1));
		checkDovetail(false, 4, RP(0, 1, 2, 1)); // not overlapping but within margin
	}
	@Test
	public void isDovetailing_should_allow_margin() {
		checkDovetail(false, 1, RP(0, 10, 8, 4));
		checkDovetail(true, 1, RP(0, 10, 9, 4));
		checkDovetail(true, 1, RP(0, 10, 10, 4));
		checkDovetail(true, 1, RP(0, 10, 11, 4));
		checkDovetail(false, 1, RP(0, 10, 12, 4));
	}
	@Test
	public void isDovetailing_should_require_both_mapped() {
		SAMRecord[] r = RP(0, 1, 1, 1);
		r[0].setReadUnmappedFlag(true);
		r[1].setMateUnmappedFlag(true);
		checkDovetail(false, 4, r);
		
		r = RP(0, 1, 1, 1);
		r[0].setMateUnmappedFlag(true);
		r[1].setReadUnmappedFlag(true);
		checkDovetail(false, 4, r);
	}
	@Test
	public void isDovetailing_should_consider_clipping_direction() {
		//           >>>>>>>>>>SSSSSSSSSS   ----> read 1
		// SSSSSSSSSS<<<<<<<<<<             <---- read 2
		//
		// classic dovetail with untrimmed adapters
		SAMRecord[] rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10M10S");
		rp[1].setCigarString("10S10M");
		clean(rp[0], rp[1]);
		assertTrue(SAMRecordUtil.isDovetailing(rp[0], PairOrientation.FR, 0));
		assertTrue(SAMRecordUtil.isDovetailing(rp[1], PairOrientation.FR, 0));
				
		
		// >>>>>>>>>>SSSSSSSSSS   ----> read 1
		// <<<<<<<<<<SSSSSSSSSS   <---- read 2
		// looks like a dovetail, but is not
		// since the read2 soft clip is on the wrong end of the read
		rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10S10M"); // <-- definitely keep this one
		rp[1].setCigarString("10S10M"); // ideally keep this one too, but we need to know about the mate cigar for that
		clean(rp[0], rp[1]);
		assertFalse(SAMRecordUtil.isDovetailing(rp[0], PairOrientation.FR, 0));
		assertFalse(SAMRecordUtil.isDovetailing(rp[1], PairOrientation.FR, 0));
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
	public void entropy_of_unspecified_sequence_should_be_4() {
		SAMRecord r = Read(0, 1, "3S1M");
		r.setReadBases(null);
		assertEquals(4, SAMRecordUtil.entropy(r), 0);
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
		read3.setReadBases(B("CGGTT"));
		read3.setBaseQualityString("54321");
		SAMRecordUtil.softenHardClips(ImmutableList.of(read, read3));
		assertEquals("AACCGACGTA", read.getReadString());
		assertEquals("1234567890", read.getBaseQualityString());
		assertEquals("TACGTCGGTT", read3.getReadString());
		assertEquals("0987654321", read3.getBaseQualityString());
	}
	
	@Test
	public void softenHardClips_should_consider_strand_negative_first() {
		SAMRecord read =  Read(0, 1, "5M5S");
		read.setReadBases(B("AACCGACGTA"));
		read.setBaseQualityString("1234567890");
		SAMRecord read3 =  Read(0, 1, "5H5M");
		read3.setReadNegativeStrandFlag(true);
		read3.setReadBases(B("CGGTT"));
		read3.setBaseQualityString("54321");
		SAMRecordUtil.softenHardClips(ImmutableList.of(read3, read));
		assertEquals("AACCGACGTA", read.getReadString());
		assertEquals("1234567890", read.getBaseQualityString());
		assertEquals("TACGTCGGTT", read3.getReadString());
		assertEquals("0987654321", read3.getBaseQualityString());
	}
	
	@Test
	public void softenHardClips_should_merge_sequences() {
		SAMRecord read =  Read(0, 1, "3M6H");
		read.setReadBases(B("ACT"));
		read.setBaseQualityString("123");
		SAMRecord read2 =  Read(0, 1, "3H3M3H");
		read2.setReadBases(B("CCA"));
		read2.setBaseQualityString("456");
		SAMRecord read3 =  Read(0, 1, "6H3M");
		read3.setReadBases(B("TGT"));
		read3.setBaseQualityString("789");
		SAMRecordUtil.softenHardClips(ImmutableList.of(read, read2, read3));
		assertEquals("ACTCCATGT", read.getReadString());
		assertEquals("ACTCCATGT", read2.getReadString());
		assertEquals("ACTCCATGT", read3.getReadString());
		assertEquals("123456789", read.getBaseQualityString());
		assertEquals("123456789", read2.getBaseQualityString());
		assertEquals("123456789", read3.getBaseQualityString());
		assertEquals("3M6S", read.getCigarString());
		assertEquals("3S3M3S", read2.getCigarString());
		assertEquals("6S3M", read3.getCigarString());
	}
	
	@Test
	public void softenHardClips_should_N_pad_if_base_cannot_be_determined() {
		SAMRecord read =  Read(0, 1, "5H5M5S");
		read.setReadBases(B("AACCGACGTA"));
		read.setBaseQualityString("1234567890");
		SAMRecord read3 =  Read(0, 1, "10H5M");
		read3.setReadBases(B("ACGTA"));
		read3.setBaseQualityString("67890");
		SAMRecordUtil.softenHardClips(ImmutableList.of(read, read3));
		assertEquals("NNNNNAACCGACGTA", read.getReadString());
		assertEquals("NNNNNAACCGACGTA", read3.getReadString());
		assertEquals("!!!!!1234567890", read.getBaseQualityString());
		assertEquals("!!!!!1234567890", read3.getBaseQualityString());
		assertEquals("10S5M", read3.getCigarString());
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
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1, read5), ImmutableSet.of(SAMTag.FI.name(), SAMTag.TC.name()), false, false, false, false, false,false);
		
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
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1, read5), ImmutableSet.of(SAMTag.R2.name(), SAMTag.Q2.name()), false, false, false, false, false,false);
		
		assertEquals("AGT", read0.getStringAttribute("R2"));
		assertEquals(null, read1.getStringAttribute("R2"));
		assertEquals("AAC", read5.getStringAttribute("R2"));
		
		assertEquals("654", read0.getStringAttribute("Q2"));
		assertEquals(null, read1.getStringAttribute("Q2"));
		assertEquals("123", read5.getStringAttribute("Q2"));
	}
	@Test
	public void calculateTemplateTags_R2_Q2_should_not_write_self() {
		SAMRecord read0 = Read(0, 1, "3M");
		read0.setReadPairedFlag(true);
		read0.setFirstOfPairFlag(true);
		read0.setReadBases(B("AAC"));
		read0.setBaseQualityString("123");
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0), ImmutableSet.of(SAMTag.R2.name(), SAMTag.Q2.name()), false, false, false, false, false,false);
		
		assertNull(read0.getAttribute("R2"));
		assertNull(read0.getAttribute("Q2"));
	}
	@Test
	public void calculateTemplateTags_should_handle_unpaired_reads() {
		SAMRecord read = Read(0, 1, "3M");
		SAMRecord anotherRead = Read(0, 1, "3M");
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read, anotherRead), SAMRecordUtil.TEMPLATE_TAGS, true, true, false, false, false,true);
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
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(SAMTag.R2.name(), SAMTag.Q2.name()), false, false, false, false, false,false);
		
		assertEquals("AGT", read0.getStringAttribute("R2"));
		assertEquals("AAC", read1.getStringAttribute("R2"));
		
		assertEquals("654", read0.getStringAttribute("Q2"));
		assertEquals("123", read1.getStringAttribute("Q2"));
	}
	@Test
	public void calculateTemplateTags_SA_should_ignore_missing_NM() {
		SAMRecord read0 = withMapq(10, withSequence("ANNNNN", Read(0, 1, "1M5S")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		
		read0.setAttribute("NM", null);
		read1.setAttribute("NM", null);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(SAMTag.SA.name()), false, false, false, true, false,false);
		
		assertEquals("polyA,2,+,1H1M4H,11,", read0.getStringAttribute("SA"));
		assertEquals("polyA,1,+,1M5S,10,", read1.getStringAttribute("SA"));
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
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(SAMTag.SA.name()), false, false, false, false, false,false);
		
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
		
		read0.setSecondaryAlignment(true);
		read1.setSecondaryAlignment(true);
		read2.setSecondaryAlignment(false);
		read3.setSecondaryAlignment(true);
		
		read0.setAttribute("NM", 0);
		read1.setAttribute("NM", 1);
		read2.setAttribute("NM", 2);
		read3.setAttribute("NM", 3);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1, read2, read3, alt), ImmutableSet.of(SAMTag.SA.name()), false, false, false, true, false,false);
		
		assertEquals(null, alt.getStringAttribute("SA"));
		// "0,4,+,3H1M2H,12,2,;0,1,+,1M5H,10,0,;0,2,+,1H1M4H,11,1,;0,5,+,4S2M,13,3,"
		//   read2                  read0            read1          read3
		assertEquals("polyA,4,+,3H1M2H,12,2;polyA,2,+,1H1M4H,11,1;polyA,5,+,4S2M,13,3", read0.getStringAttribute("SA"));
		assertEquals("polyA,4,+,3H1M2H,12,2;polyA,1,+,1M5H,10,0;polyA,5,+,4S2M,13,3", read1.getStringAttribute("SA"));
		assertEquals("polyA,1,+,1M5H,10,0;polyA,2,+,1H1M4H,11,1;polyA,5,+,4S2M,13,3", read2.getStringAttribute("SA"));
		assertEquals("polyA,4,+,3H1M2H,12,2;polyA,1,+,1M5H,10,0;polyA,2,+,1H1M4H,11,1", read3.getStringAttribute("SA"));
	}
	@Test
	public void calculateMultimappingTags_should_use_circular_pointers() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(1, 2, "1H1M4H")))[0];
		SAMRecord read2 = withMapq(12, withSequence("G", Read(3, 4, "3H1M2H")))[0];
		SAMRecord read3 = withMapq(13, withSequence("ACGTAA", Read(2, 5, "4S2M")))[0];
		
		read0.setSupplementaryAlignmentFlag(true);
		read1.setSupplementaryAlignmentFlag(true);
		read2.setSupplementaryAlignmentFlag(false); // read 2 is the 'primary' chimeric alignment record
		read3.setSupplementaryAlignmentFlag(true);
		
		SAMRecordUtil.calculateMultimappingTags(ImmutableSet.of(SAMTag.IH.name(), SAMTag.HI.name(), SAMTag.CC.name(), SAMTag.CP.name()),
				Lists.newArrayList(ImmutableList.of(read0, read1, read2, read3)));
		
		assertEquals(2, (int)read0.getIntegerAttribute("CP"));
		assertEquals(5, (int)read1.getIntegerAttribute("CP"));
		assertEquals(1, (int)read2.getIntegerAttribute("CP"));
		assertEquals(4, (int)read3.getIntegerAttribute("CP"));
				
		assertEquals("polyACGT", read0.getStringAttribute("CC"));
		assertEquals("random", read1.getStringAttribute("CC"));
		assertEquals("polyA", read2.getStringAttribute("CC"));
		assertEquals("Npower2", read3.getStringAttribute("CC"));
		
	}
	@Test
	public void calculateMultimappingTags_HI_should_number_reads_by_position() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(1, 2, "1H1M4H")))[0];
		SAMRecord read2 = withMapq(12, withSequence("G", Read(3, 4, "3H1M2H")))[0];
		SAMRecord read3 = withMapq(13, withSequence("ACGTAA", Read(2, 5, "4S2M")))[0];
		
		read0.setSupplementaryAlignmentFlag(true);
		read1.setSupplementaryAlignmentFlag(true);
		read2.setSupplementaryAlignmentFlag(false); // read 2 is the 'primary' chimeric alignment record
		read3.setSupplementaryAlignmentFlag(true);
		
		SAMRecordUtil.calculateMultimappingTags(ImmutableSet.of(SAMTag.IH.name(), SAMTag.HI.name(), SAMTag.CC.name(), SAMTag.CP.name()),
				Lists.newArrayList(ImmutableList.of(read0, read1, read2, read3)));
		
		assertEquals(0, (int)read0.getIntegerAttribute("HI"));
		assertEquals(1, (int)read1.getIntegerAttribute("HI"));
		assertEquals(3, (int)read2.getIntegerAttribute("HI"));
		assertEquals(2, (int)read3.getIntegerAttribute("HI"));
	}
	@Test
	public void calculateMultimappingTags_IH_count_reads_TODO_what_to_do_about_supplimentary_alignments() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(1, 2, "1H1M4H")))[0];
		SAMRecord read2 = withMapq(12, withSequence("G", Read(3, 4, "3H1M2H")))[0];
		SAMRecord read3 = withMapq(13, withSequence("ACGTAA", Read(2, 5, "4S2M")))[0];
		
		read0.setSupplementaryAlignmentFlag(true);
		read1.setSupplementaryAlignmentFlag(false);
		read2.setSupplementaryAlignmentFlag(false); // read 2 is the 'primary' chimeric alignment record
		read3.setSupplementaryAlignmentFlag(true);
		
		SAMRecordUtil.calculateMultimappingTags(ImmutableSet.of(SAMTag.IH.name(), SAMTag.HI.name(), SAMTag.CC.name(), SAMTag.CP.name()),
				Lists.newArrayList(ImmutableList.of(read0, read1, read2, read3)));
		
		assertEquals(4, (int)read0.getIntegerAttribute("IH"));
		assertEquals(4, (int)read1.getIntegerAttribute("IH"));
		assertEquals(4, (int)read2.getIntegerAttribute("IH"));
		assertEquals(4, (int)read3.getIntegerAttribute("IH"));
	}
	@Test
	public void lowMapqToUnmapped_should_convert_below_threshold() {
		assertTrue(SAMRecordUtil.lowMapqToUnmapped(withMapq(0, Read(0, 3, "4S5M"))[0], 2).getReadUnmappedFlag());
		assertTrue(SAMRecordUtil.lowMapqToUnmapped(withMapq(1, Read(0, 3, "4S5M"))[0], 2).getReadUnmappedFlag());
		assertFalse(SAMRecordUtil.lowMapqToUnmapped(withMapq(2, Read(0, 3, "4S5M"))[0], 2).getReadUnmappedFlag());
		assertFalse(SAMRecordUtil.lowMapqToUnmapped(withMapq(3, Read(0, 3, "4S5M"))[0], 2).getReadUnmappedFlag());
		assertFalse(SAMRecordUtil.lowMapqToUnmapped(withMapq(255, Read(0, 3, "4S5M"))[0], 2).getReadUnmappedFlag());
	}
	@Test
	public void lowMapqToUnmapped_should_adjust_chimeric_fragments() {
		assertEquals("polyA,2,+,5M,2,0", SAMRecordUtil.lowMapqToUnmapped(withAttr("SA", "polyA,1,+,5M,1,0;polyA,2,+,5M,2,0", Read(0, 3, "4S5M"))[0], 2).getStringAttribute("SA"));
	}
	/**
	 * Since the scoring is based on the length of the primary read soft clip,
	 * Unmapping the primary without also unmapping all supplementary alignments
	 * will cause asymmetry in the scoring.
	 */
	@Test
	public void lowMapqToUnmapped_should_consider_all_fragments_unmapped_if_primary_unmapped() {
		SAMRecord r = withAttr("SA", "polyA,1,+,5M,1,0;polyA,2,+,5M,2,0", Read(0, 3, "4S5M"))[0];
		r.setSupplementaryAlignmentFlag(false);
		assertFalse(SAMRecordUtil.lowMapqToUnmapped(r, 2).getReadUnmappedFlag());
		
		r = withAttr("SA", "polyA,1,+,5M,1,0;polyA,2,+,5M,2,0", Read(0, 3, "4S5M"))[0];
		r.setSupplementaryAlignmentFlag(true);
		assertTrue(SAMRecordUtil.lowMapqToUnmapped(r, 2).getReadUnmappedFlag());
	}
	@Test
	public void lowMapqToUnmapped_should_adjust_mate() {
		assertTrue(SAMRecordUtil.lowMapqToUnmapped(withAttr("MQ", 0, RP(0, 1, 1))[0], 2).getMateUnmappedFlag());
		assertTrue(SAMRecordUtil.lowMapqToUnmapped(withAttr("MQ", 1, RP(0, 1, 1))[0], 2).getMateUnmappedFlag());
		assertFalse(SAMRecordUtil.lowMapqToUnmapped(withAttr("MQ", 2, RP(0, 1, 1))[0], 2).getMateUnmappedFlag());
		assertFalse(SAMRecordUtil.lowMapqToUnmapped(withAttr("MQ", 3, RP(0, 1, 1))[0], 2).getMateUnmappedFlag());
	}
	@Test
	public void should_fix_mate_information() {
		SAMRecord r1 = Read(0, 1, "10M1I1M");
		SAMRecord r2 = Read(2, 3, "10M5S");
		r2.setAttribute("FI", 1);
		r1.setReadName("r");
		r2.setReadName("r");
		r1.setMappingQuality(5);
		r2.setMappingQuality(6);
		r2.setReadNegativeStrandFlag(true);
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(SAMTag.MC.name(), SAMTag.MQ.name()), false, true, false, false, false, false);
		
		assertEquals(2, (int)r1.getMateReferenceIndex());
		assertEquals(3, r1.getMateAlignmentStart());
		assertEquals(0, (int)r2.getMateReferenceIndex());
		assertEquals(1, r2.getMateAlignmentStart());
		assertTrue(r1.getMateNegativeStrandFlag());
		assertFalse(r2.getMateNegativeStrandFlag());
		assertFalse(r1.getMateUnmappedFlag());
		assertFalse(r2.getMateUnmappedFlag());
		assertEquals("10M5S", r1.getStringAttribute("MC"));
		assertEquals("10M1I1M", r2.getStringAttribute("MC"));
		assertEquals(6, (int)r1.getIntegerAttribute("MQ"));
		assertEquals(5, (int)r2.getIntegerAttribute("MQ"));
	}
	@Test
	public void fix_mate_information_should_leave_unpaired() {
		SAMRecord r = Read(0, 15, "5M6S");
		assertFalse(r.getReadPairedFlag());
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r), ImmutableSet.of(SAMTag.MC.name(), SAMTag.MQ.name()), false, true, false, false, false, false);
		assertFalse(r.getReadPairedFlag());
	}@Test
	public void fix_duplicate_should_mark_all_reads_if_any_are_duplicates() {
		SAMRecord[] rp = RP(0, 1, 100);
		rp[1].setDuplicateReadFlag(true);
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(rp[0], rp[1]), ImmutableSet.of(), false, false, true, false, false, false);
		Assert.assertTrue(Arrays.stream(rp).allMatch(r -> r.getDuplicateReadFlag()));
	}
	@Test
	public void recalculateSupplementary_should_convert_to_secondary_SA_to_supplementary() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		
		read0.setAttribute("SA", "polyA,2,+,1H1M4H,11,0");
		read1.setAttribute("SA", "polyA,1,+,1M5H,10,0");
		read0.setSecondaryAlignment(true);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(), false, false, false, false, false,true);
		
		Assert.assertTrue(read0.getSupplementaryAlignmentFlag());
		Assert.assertFalse(read1.getSupplementaryAlignmentFlag());
	}
	@Test
	public void recalculateSupplementary_should_force_primary_to_supp() {
		SAMRecord read0 = withMapq(10, withSequence("ANNNNN", Read(0, 1, "1M5S")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		
		read0.setAttribute("SA", "polyA,2,+,1H1M4H,11,0");
		read1.setAttribute("SA", "polyA,1,+,1M5S,10,0");
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(), false, false, false, false, false,true);
		
		Assert.assertFalse(read0.getSupplementaryAlignmentFlag());
		Assert.assertTrue(read1.getSupplementaryAlignmentFlag());
	}
	@Test
	public void recalculateSupplementary_should_not_convert_unmapped_read_to_supplementary() {
		SAMRecord read0 = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1 = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		
		read0.setAttribute("SA", "polyA,2,+,1H1M4H,11,0");
		read1.setAttribute("SA", "polyA,1,+,1M5H,10,0");
		
		read0.setReadUnmappedFlag(true);
		read1.setReadUnmappedFlag(true);
		
		read1.setSecondaryAlignment(true);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read0, read1), ImmutableSet.of(), false, false, false, false, false,true);
		
		Assert.assertFalse(read0.getSupplementaryAlignmentFlag());
		Assert.assertFalse(read1.getSupplementaryAlignmentFlag());
	}
	@Test
	public void recalculateSupplementary_should_recalculate_supplementary_based_on_SA_tag() {
		// primary split alignment
		SAMRecord read1a = withMapq(10, withSequence("A", Read(0, 1, "1M5H")))[0];
		SAMRecord read1b = withMapq(11, withSequence("C", Read(0, 2, "1H1M4H")))[0];
		
		// secondary
		SAMRecord read2 = withMapq(11, withSequence("CTG", Read(0, 2, "3M")))[0];
		
		// secondary split alignment
		SAMRecord read3a = withMapq(10, withSequence("AA", Read(1, 1, "2M6H")))[0];
		SAMRecord read3b = withMapq(11, withSequence("CGTCGT", Read(1, 2, "2H3M3S")))[0];
		//SAMRecord read3c = withMapq(11, withSequence("C", Read(1, 2, "2H1M4H")))[0];
		
		read1a.setAttribute("SA", "polyA,2,+,1H1M4H,11,0");
		read1b.setAttribute("SA", "polyA,1,+,1M5H,10,0");
		
		read3a.setAttribute("SA", "polyACGT,2,+,2H3M3S,11,0");
		read3b.setAttribute("SA", "polyACGT,1,+,2M6H,10,0");
		
		read1a.setSupplementaryAlignmentFlag(false);
		read1b.setSupplementaryAlignmentFlag(true);
		read2.setSupplementaryAlignmentFlag(false);
		read3a.setSupplementaryAlignmentFlag(false);
		read3b.setSupplementaryAlignmentFlag(false);
		
		read1a.setSecondaryAlignment(false);
		read1b.setSecondaryAlignment(false);
		read2.setSecondaryAlignment(true);
		read3a.setSecondaryAlignment(true);
		read3b.setSecondaryAlignment(true);
		
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(read1a, read1b, read2, read3a, read3b), ImmutableSet.of(), false, false, false, false, false,true);
		
		Assert.assertFalse(read1a.getSupplementaryAlignmentFlag());
		Assert.assertTrue(read1b.getSupplementaryAlignmentFlag());
		Assert.assertFalse(read2.getSupplementaryAlignmentFlag());
		Assert.assertTrue(read3a.getSupplementaryAlignmentFlag());
		Assert.assertFalse(read3b.getSupplementaryAlignmentFlag());
		
		Assert.assertFalse(read1a.isSecondaryAlignment());
		Assert.assertFalse(read1b.isSecondaryAlignment());
		Assert.assertTrue(read2.isSecondaryAlignment());
		Assert.assertTrue(read3a.isSecondaryAlignment());
		Assert.assertTrue(read3b.isSecondaryAlignment());
	}
	@Test
	public void unclipExactReferenceMatches_should_adjust_matches() {
		SAMRecord r = Read(0, 10, "5S5M5S");
		SAMRecordUtil.unclipExactReferenceMatches(SMALL_FA, r);
		assertEquals(5, r.getAlignmentStart());
		assertEquals("15M", r.getCigarString());
		// 12345678
		// ACGTACGT
		// gcgtacgg
		// SSSMSSSS
		// SMMMMMMS
		r = withSequence("gcgtacgg", Read(1, 4, "3S1M4S"))[0];
		SAMRecordUtil.unclipExactReferenceMatches(SMALL_FA, r);
		assertEquals(2, r.getAlignmentStart());
		assertEquals("1S6M1S", r.getCigarString());
	}
	@Test
	public void unclipExactReferenceMatches_should_allow_hard_clipping() {
		SAMRecord r = Read(0, 10, "10H5S5M5S15H");
		SAMRecordUtil.unclipExactReferenceMatches(SMALL_FA, r);
		assertEquals(5, r.getAlignmentStart());
		assertEquals("10H15M15H", r.getCigarString());
	}
	@Test
	public void unclipExactReferenceMatches_should_not_overrun_reference() {
		SAMRecord r = withSequence("NN", Read(0, 1, "1S1M"))[0];
		SAMRecordUtil.unclipExactReferenceMatches(SMALL_FA, r);
		assertEquals(1, r.getAlignmentStart());
	}
	@Test
	public void unclipExactReferenceMatches_should_unclip_entire_read() {
		String refStr = "TAAATTGGAACACTATACCAAAACATTAACCAGCATAGCAGTATATAAGGTTAAACATTAAATAACCCCTGGCTTAACTAACTCTCCAATTGCACTTTCTATAAGTAATTGTTGTTTAGACTTTATTAATTCAGATGTTTCAGACATGTCTTATATACACAAGAGAATTTCATTTCTCTTT";
		String readStr = "AAATTGGAACACTATACCAAAACATTAACCAGCATAGCAGTATATAAGGTTAAACATTAAATAACCCCTGGCTTAACTAACTCTCCAATTGCACTTTCTATAAGTAATTGTTGTTTAGACTTTATTAATTC";
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "Contig" }, new byte[][] { B(refStr) });
		SAMRecord r = new SAMRecord(new SAMFileHeader());
		r.getHeader().setSequenceDictionary(ref.getSequenceDictionary());
		r.setReferenceIndex(0);
		r.setCigarString("97M34S");
		r.setAlignmentStart(2);
		r.setReadNegativeStrandFlag(false);
		r.setReadBases(B(readStr));
		SAMRecordUtil.unclipExactReferenceMatches(ref, r);
		assertEquals("131M", r.getCigarString());
	}
	@Test
	public void hardClipToN_should_pad_with_zero_qual_Ns() {
		SAMRecord r = withQual(new byte[] {1, 2}, withSequence("GT", Read(0, 1, "1H2M3H")))[0];
		SAMRecordUtil.hardClipToN(r);
		assertEquals("1S2M3S", r.getCigarString());
		assertEquals("NGTNNN", S(r.getReadBases()));
		assertArrayEquals(new byte[] {0, 1, 2, 0,0,0}, r.getBaseQualities());
	}
	@Test
	public void getEffectiveMapq_should_use_mapq() {
		assertEquals(0, SAMRecordUtil.getEffectiveMapq(withMapq(0, Read(0, 1, "1M"))[0], 1), 0);
		assertEquals(1, SAMRecordUtil.getEffectiveMapq(withMapq(1, Read(0, 1, "1M"))[0], 2), 0);
	}
	@Test
	public void getEffectiveMapq_should_be_0_for_unmapped_reads() {
		SAMRecord r = withMapq(1, Read(0, 1, "1M"))[0];
		r.setReadUnmappedFlag(true);
		assertEquals(0, SAMRecordUtil.getEffectiveMapq(r, 1), 0);
	}
	@Test
	public void getEffectiveMapq_should_use_alternate_alignment_location_count() {
		assertEquals(3.010299957, SAMRecordUtil.getEffectiveMapq(withAttr("IH", 2, withMapq(10, Read(0, 1, "1M")))[0], 1), 0.00001);
		assertEquals(0.457574905, SAMRecordUtil.getEffectiveMapq(withAttr("NH", 10, withMapq(15, Read(0, 1, "1M")))[0], 1), 0.00001);
	}
	@Test
	public void getEffectiveMapq_should_use_fallback() {
		assertEquals(1.5, SAMRecordUtil.getEffectiveMapq(withMapq(SAMRecord.UNKNOWN_MAPPING_QUALITY, Read(0, 1, "1M"))[0], 1.5), 0);
	}
	@Test
	public void restoreHardClips_should_not_adjust_alignment() {
		SAMRecord r1 = Read(0, 100, "48M52S");
		SAMRecord r2 = Read(0, 100, "32H48M20H");
		SAMRecord r3 = Read(0, 100, "64H36M");
		r1.setAttribute("SA", new ChimericAlignment(r2).toString() + ";" + new ChimericAlignment(r3).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString() + ";" + new ChimericAlignment(r3).toString());
		r3.setAttribute("SA", new ChimericAlignment(r1).toString() + ";" + new ChimericAlignment(r2).toString());
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2, r3), ImmutableSet.of(), true, false, false, false, false,false);
		Assert.assertEquals("32S48M20S", r2.getCigarString());
		Assert.assertEquals("64S36M", r3.getCigarString());
	}
	@Test
	public void should_return_overlapping_bases() {
		Assert.assertEquals(0, SAMRecordUtil.overlappingBases(Read(0, 1, "10M"), Read(0, 20, "10M")));
		Assert.assertEquals(0, SAMRecordUtil.overlappingBases(Read(0, 1, "10M"), onNegative(Read(1, 1, "10M"))[0]));
		Assert.assertEquals(10, SAMRecordUtil.overlappingBases(Read(0, 1, "5M5D10M"), Read(0, 1, "15M")));
	}
	@Test
	public void should_not_overlap_on_different_contigs() {
		Assert.assertEquals(0, SAMRecordUtil.overlappingBases(Read(0, 1, "10M"), Read(1, 1, "10M")));
	}
	@Test
	public void should_not_overlap_on_opposite_strand() {
		Assert.assertEquals(0, SAMRecordUtil.overlappingBases(Read(0, 1, "10M"), onNegative(Read(1, 1, "10M"))[0]));
	}
	@Test
	public void getReadLengthIncludingHardClipping() {
		Assert.assertEquals(6, SAMRecordUtil.getReadLengthIncludingHardClipping(Read(0, 1, "1H2M3H")));
	}
	@Test
	public void adjustAlignmentBounds_should_adjust_bounds() {
		Assert.assertEquals("1S7M2S", SAMRecordUtil.adjustAlignmentBounds(Read(0, 1, "10M"), -1, -2).getCigarString());
		Assert.assertEquals("1H8S2M", SAMRecordUtil.adjustAlignmentBounds(Read(0, 1, "1H2S3M4D5M"), -6, 0).getCigarString());
		Assert.assertEquals("2M8S1H", SAMRecordUtil.adjustAlignmentBounds(Read(0, 1, "5M4D3M2S1H"), 0, -6).getCigarString());
	}
	@Test
	public void adjustAlignmentBounds_should_convert_negative_to_soft_clip() {
		Assert.assertEquals("1X1=1X1=1X5S", SAMRecordUtil.adjustAlignmentBounds(Read(0, 1, "1X1=1X1=1X1=1X1=1X1="), 0, -5).getCigarString());
	}
	@Test
	public void getEndClipLength_should_consider_end_clipped_if_no_reads_mapped() {
		Assert.assertEquals(5, SAMRecordUtil.getEndClipLength(Read(0, 1, "5S")));
	}
	@Test
	public void fixSA_should_regenerate_SA_tag() {
		SAMRecord r1 = Read(0, 100, "10M20S");
		SAMRecord r2 = Read(0, 200, "10S10M10S");
		SAMRecord r3 = Read(0, 300, "20S10M");
		r1.setAttribute("SA", new ChimericAlignment(r1).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString());
		r3.setAttribute("SA", new ChimericAlignment(r1).toString());
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2, r3), ImmutableSet.of(SAMTag.SA.name()), false, false, false, true, false,false);
		Assert.assertEquals("polyA,200,+,10S10M10S,10,0;polyA,300,+,20S10M,10,0", r1.getStringAttribute("SA"));
		Assert.assertEquals("polyA,100,+,10M20S,10,0;polyA,300,+,20S10M,10,0", r2.getStringAttribute("SA"));
		Assert.assertEquals("polyA,100,+,10M20S,10,0;polyA,200,+,10S10M10S,10,0", r3.getStringAttribute("SA"));
	}
	@Test
	public void fixSA_should_consider_first_primary_record_as_primary() {
		SAMRecord r1 = Read(0, 100, "10M20S");
		SAMRecord r2 = Read(0, 200, "10S10M10S");
		SAMRecord r3 = Read(0, 300, "20S10M");
		r1.setAttribute("SA", new ChimericAlignment(r1).toString());
		r2.setAttribute("SA", new ChimericAlignment(r1).toString());
		r3.setAttribute("SA", new ChimericAlignment(r1).toString());
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2, r3), ImmutableSet.of(SAMTag.SA.name()), false, false, false, true, false,true);
		Assert.assertFalse(r1.getSupplementaryAlignmentFlag());
		Assert.assertFalse(r1.isSecondaryAlignment());
		Assert.assertTrue(r2.getSupplementaryAlignmentFlag());
		Assert.assertFalse(r2.isSecondaryAlignment());
		Assert.assertTrue(r3.getSupplementaryAlignmentFlag());
		Assert.assertFalse(r3.isSecondaryAlignment());
	}
	@Test
	public void fixSA_convert_overlapping_read_to_split_read_alignment() {
		SAMRecord r1 = Read(0, 100, "70M30S");
		SAMRecord r2 = Read(0, 200, "30S70M");
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(SAMTag.SA.name()), false, false, false, true, false,true);
		Assert.assertFalse(r1.getSupplementaryAlignmentFlag());
		Assert.assertFalse(r1.isSecondaryAlignment());
		Assert.assertTrue(r2.getSupplementaryAlignmentFlag());
		Assert.assertFalse(r2.isSecondaryAlignment());
		Assert.assertEquals("polyA,200,+,30S70M,10,0", r1.getStringAttribute("SA"));
		Assert.assertEquals("polyA,100,+,70M30S,10,0", r2.getStringAttribute("SA"));
	}
	@Test
	public void fixSA_should_remove_duplicate_alignments_starting_at_same_read_offset() {
		SAMRecord r1 = Read(0, 100, "70M30S");
		SAMRecord r2 = Read(0, 200, "50M50S");
		ArrayList list = Lists.newArrayList(r1, r2);
		SAMRecordUtil.calculateTemplateTags(list, ImmutableSet.of(SAMTag.SA.name()), false, false, false, true, false,true);
	}
	@Test
	public void fixTruncated_should_add_hard_clipping_to_shorter_read_to_match_longer() {
		// ACGTTTGGAAT
		//    TTTGG
		// HHHsmmmmHHH
		SAMRecord r1 = withSequence("ACGTTTGGAAT", Read(0, 100, "11M"))[0];
		SAMRecord r2 = withSequence("TTTGG", Read(0, 100, "1S4M1H"))[0];
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(), false, false, false, false, true,false);
		Assert.assertEquals("3H1S4M3H", r2.getCigarString());
	}
	@Test
	public void fixTruncated_should_consider_alignment_strand() {
		//  TTCCAAACGT
		//    CCAAA
		//  HHsmmmmHHH
		SAMRecord r1 = withSequence("ACGTTTGGAA", Read(0, 100, "10M"))[0];
		SAMRecord r2 = onNegative(withSequence("CCAAA", Read(0, 100, "1S4M1H")))[0];
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(), false, false, false, false, true,false);
		Assert.assertEquals("2H1S4M3H", r2.getCigarString());
	}
	@Test
	public void fixTruncated_should_match_at_upper_bound() {
		SAMRecord r1 = withSequence("ACGTTTGGAA", Read(0, 100, "10M"))[0];
		SAMRecord r2 = withSequence("TGGAA", Read(0, 100, "1S4M"))[0];
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(), false, false, false, false, true,false);
		Assert.assertEquals("5H1S4M", r2.getCigarString());
	}
	@Test
	public void fixTruncated_should_match_at_lower_bound() {
		SAMRecord r1 = withSequence("ACGTTTGGAA", Read(0, 100, "10M"))[0];
		SAMRecord r2 = withSequence("ACGTT", Read(0, 100, "1S4M"))[0];
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(), false, false, false, false, true,false);
		Assert.assertEquals("1S4M5H", r2.getCigarString());
	}
	@Test
	public void fixTruncated_should_not_update_if_sequence_does_not_match() {
		SAMRecord r1 = withSequence("AGGTTTGGAA", Read(0, 100, "10M"))[0];
		SAMRecord r2 = withSequence("ACGTT", Read(0, 100, "1S4M1H"))[0];
		SAMRecordUtil.calculateTemplateTags(ImmutableList.of(r1, r2), ImmutableSet.of(), false, false, false, false, true,false);
		Assert.assertEquals("1S4M1H", r2.getCigarString());
	}
	@Test
	public void forceValidContigBounds_should_force_to_contig_start() {
		// 2 1
		// ^ ^012345
		// MIMDDDMMM
		SAMRecord r = Read(0, -2, "1M1I1M3D3M");
		SAMRecordUtil.forceValidContigBounds(r, SMALL_FA.getSequenceDictionary());
		assertEquals(3, r.getAlignmentStart());
		assertEquals("3S3M", r.getCigarString());


		r = Read(0, -1, "10M");
		SAMRecordUtil.forceValidContigBounds(r, SMALL_FA.getSequenceDictionary());
		assertEquals(1, r.getAlignmentStart());
		assertEquals("2S8M", r.getCigarString());

		r = Read(0, -1, "2I10M");
		SAMRecordUtil.forceValidContigBounds(r, SMALL_FA.getSequenceDictionary());
		assertEquals(1, r.getAlignmentStart());
		assertEquals("4S8M", r.getCigarString());
	}
	@Test
	public void forceValidContigBounds_should_force_to_contig_end() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] {"contig"}, new byte[][] { new byte[] {0,1,2,3,4,5,6,7}});
		//         90123456
		// 12345678
		//   MMMDDDDIMDM
		//   MMM    SS S
		SAMRecord r = Read(0, 3, "3M4D1I1M1D1M");
		SAMRecordUtil.forceValidContigBounds(r, ref.getSequenceDictionary());
		assertEquals(5, r.getAlignmentEnd());
		assertEquals("3M3S", r.getCigarString());
	}
	//@Test
	@Ignore("Handled downstream instead so the dovetailing check doesn't need to chase split read alignments")
	public void isDovetailing_should_also_test_primary_alignment() throws CloneNotSupportedException {
		SAMRecord[] rp = RP(2, 100, 100, 100);
		rp[0].setCigarString("75M25S");
		rp[1].setCigarString("25S75M");
		rp[0].setAttribute("MC", rp[1].getCigarString());
		rp[1].setAttribute("MC", rp[0].getCigarString());
		SAMRecord supp = (SAMRecord) rp[0].clone();
		supp.setSupplementaryAlignmentFlag(true);
		supp.setAlignmentStart(200);
		supp.setCigarString("75S25M");
		supp.setAttribute("SA", new ChimericAlignment(rp[0]).toString());
		rp[0].setAttribute("SA", new ChimericAlignment(supp).toString());
		Assert.assertTrue(SAMRecordUtil.isDovetailing(rp[0], PairOrientation.FR, 0));
		Assert.assertTrue(SAMRecordUtil.isDovetailing(supp, PairOrientation.FR, 0));
	}
	@Test
	public void fixMate_should_set_to_primary() throws CloneNotSupportedException {
		SAMRecord[] rp = RP(2, 100, 100, 100);
		rp[0].setCigarString("75M25S");
		rp[1].setCigarString("25S75M");
		SAMRecord[] supp = new SAMRecord[] { (SAMRecord)rp[0].clone(), (SAMRecord)rp[1].clone()};
		supp[0].setSupplementaryAlignmentFlag(true);
		supp[1].setSupplementaryAlignmentFlag(true);
		supp[0].setAlignmentStart(200);
		supp[1].setAlignmentStart(300);
		supp[0].setCigarString("60S40M");
		supp[1].setCigarString("10M90S");

		rp[0].setMappingQuality(11);
		rp[1].setMappingQuality(12);
		supp[0].setMappingQuality(13);
		supp[1].setMappingQuality(14);
		rp[0].setAttribute("SA", new ChimericAlignment(supp[0]).toString());
		rp[1].setAttribute("SA", new ChimericAlignment(supp[1]).toString());
		supp[0].setAttribute("SA", new ChimericAlignment(rp[0]).toString());
		supp[1].setAttribute("SA", new ChimericAlignment(rp[1]).toString());

		List<List<SAMRecord>> fragment = Lists.newArrayList(
				Lists.newArrayList(supp[0], rp[0]),
				Lists.newArrayList(supp[1], rp[1]));
		SAMRecordUtil.fixMates(fragment, true, true);

		Assert.assertEquals("25S75M", rp[0].getStringAttribute("MC"));
		Assert.assertEquals("75M25S", rp[1].getStringAttribute("MC"));
		Assert.assertEquals("25S75M", supp[0].getStringAttribute("MC"));
		Assert.assertEquals("75M25S", supp[1].getStringAttribute("MC"));

		Assert.assertEquals(12, (int)rp[0].getIntegerAttribute("MQ"));
		Assert.assertEquals(11, (int)rp[1].getIntegerAttribute("MQ"));
		Assert.assertEquals(12, (int)supp[0].getIntegerAttribute("MQ"));
		Assert.assertEquals(11, (int)supp[1].getIntegerAttribute("MQ"));
	}
	@Test
	public void matchReadPairPrimaryAlignments_should_pair_primaries() throws CloneNotSupportedException {
		SAMRecord[] rp = RP(2, 100, 100, 100);
		rp[0].setCigarString("75M25S");
		rp[1].setCigarString("25S75M");
		SAMRecord[] supp = new SAMRecord[] { (SAMRecord)rp[0].clone(), (SAMRecord)rp[1].clone()};
		supp[0].setSupplementaryAlignmentFlag(true);
		supp[1].setSupplementaryAlignmentFlag(true);
		supp[0].setAlignmentStart(200);
		supp[1].setAlignmentStart(300);
		supp[0].setCigarString("60S40M");
		supp[1].setCigarString("10M90S");

		rp[0].setMappingQuality(11);
		rp[1].setMappingQuality(12);
		supp[0].setMappingQuality(13);
		supp[1].setMappingQuality(14);
		rp[0].setAttribute("SA", new ChimericAlignment(supp[0]).toString());
		rp[1].setAttribute("SA", new ChimericAlignment(supp[1]).toString());
		supp[0].setAttribute("SA", new ChimericAlignment(rp[0]).toString());
		supp[1].setAttribute("SA", new ChimericAlignment(rp[1]).toString());

		SamPairUtil.setMateInfo(rp[0], supp[1], true);
		SamPairUtil.setMateInfo(rp[1], supp[0], true);

		List<List<SAMRecord>> fragment = Lists.newArrayList(
				Lists.newArrayList(supp[0], rp[0]),
				Lists.newArrayList(supp[1], rp[1]));
		SAMRecordUtil.matchReadPairPrimaryAlignments(fragment);

		Assert.assertEquals("25S75M", rp[0].getStringAttribute("MC"));
		Assert.assertEquals("75M25S", rp[1].getStringAttribute("MC"));

		Assert.assertEquals(12, (int)rp[0].getIntegerAttribute("MQ"));
		Assert.assertEquals(11, (int)rp[1].getIntegerAttribute("MQ"));
	}
	@Test
	@Category(EColiTests.class)
	public void issue278_inconsistent_read_pair() throws FileNotFoundException {
		File ref = ReferenceTests.findReference("Escherichia_coli_bl21_de3_.ASM956v1.dna.toplevel.fa");
		SynchronousReferenceLookupAdapter reflookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
		File input = new File("src/test/resources/sanity_failure_debug/bl21_de3_.ASM956v1_S2.bam");
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(input.getParentFile(), input.getParentFile(), SAMFileWriterImpl.getDefaultMaxRecordsInRam()), ref, reflookup, null, getConfig());
		List<SAMRecord> results = Lists.newArrayList(new UngroupingIterator<>(new TemplateTagsIterator(
				TemplateTagsIterator.withGrouping(getRecords(input).iterator()),
				true, true, true, true, true, true, ImmutableSet.of(SAMTag.MC.name(), SAMTag.MQ.name()))));
		List<SAMRecord> primary = results.stream().filter(r -> !r.isSecondaryOrSupplementary()).collect(Collectors.toList());
		assertEquals(primary.get(0).getCigarString(), primary.get(1).getStringAttribute("MC"));
		assertEquals(primary.get(1).getCigarString(), primary.get(0).getStringAttribute("MC"));
	}
	@Test
	@Category(EColiTests.class)
	public void issue278_inconsistent_assembly_scoring() throws IOException {
		File ref = ReferenceTests.findReference("Escherichia_coli_bl21_de3_.ASM956v1.dna.toplevel.fa");
		SynchronousReferenceLookupAdapter reflookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
		File input = new File("src/test/resources/sanity_failure_debug/bl21_de3_.ASM956v1_assembly.bam");
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(input.getParentFile(), input.getParentFile(), SAMFileWriterImpl.getDefaultMaxRecordsInRam()), ref, reflookup, null, getConfig());
		List<SAMRecord> records = getRecords(input).stream()
				.filter(r -> r.getReadName().equals("asm0-758922"))
				.collect(Collectors.toList());
		MockSAMEvidenceSource ses = SES(pc);
		List<DirectedBreakpoint> evidence = records.stream()
				.flatMap(r -> SingleReadEvidence.createEvidence(ses, 0, r).stream())
				.filter(r -> r instanceof DirectedBreakpoint)
				.map(r -> (DirectedBreakpoint)r)
				.collect(Collectors.toList());
		Assert.assertEquals(evidence.get(0).getBreakpointQual(), evidence.get(1).getBreakpointQual(), 0);
	}
	@Test
	public void fixMate_should_not_match_to_supplementary_if_primary_exists() throws CloneNotSupportedException {
		SAMRecord[] rp = RP(0, 100, 200, 10);
		SAMRecord primary = (SAMRecord)rp[1].clone();
		rp[1].setSupplementaryAlignmentFlag(true);
		primary.setAlignmentStart(300);
		primary.getMateUnmappedFlag();
		primary.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
		primary.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);

		SAMRecordUtil.fixMates(ImmutableList.of(Lists.newArrayList(rp[0]), Lists.newArrayList(rp[1], primary)), true, true);
		assertEquals(primary.getAlignmentStart(), rp[0].getMateAlignmentStart());
	}
	@Test
	public void issue312_should_handle_multiple_overlapping_records() {
		List<SAMRecord> list = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			list.add(Read(0, i, i + "S1M" + (10-i) + "S"));
			list.add(Read(0, i, i + "S1M" + (10-i) + "S"));
		}
		list.stream().forEach(r -> {
			r.setReadName("r");
			r.setSupplementaryAlignmentFlag(true);
		});
		list.get(0).setSupplementaryAlignmentFlag(false);
		SAMRecordUtil.calculateTemplateTags(list, ImmutableSet.of(SAMTag.SA.name()), false, false,false, true, false, false);
	}
}


