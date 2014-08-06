package au.edu.wehi.idsv.sam;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import au.edu.wehi.idsv.sam.SAMRecordUtil;


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
	public void isPartOfNonReferenceReadPair() {
		assertFalse(SAMRecordUtil.isPartOfNonReferenceReadPair(Read(0, 1, "1M")));
		assertFalse(SAMRecordUtil.isPartOfNonReferenceReadPair(RP(1, 1, 1)[0]));
		assertFalse(SAMRecordUtil.isPartOfNonReferenceReadPair(RP(1, 1, 1)[1]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(OEA(1, 1, "1M", true)[0]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(OEA(1, 1, "1M", false)[0]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(OEA(1, 1, "1M", true)[1]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(OEA(1, 1, "1M", false)[1]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", true, 0, 1, "1M", true)[0]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", true, 0, 1, "1M", false)[0]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", false, 0, 1, "1M", true)[0]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", false, 0, 1, "1M", false)[0]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", true, 0, 1, "1M", true)[1]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", true, 0, 1, "1M", false)[1]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", false, 0, 1, "1M", true)[1]));
		assertTrue(SAMRecordUtil.isPartOfNonReferenceReadPair(DP(1, 1, "1M", false, 0, 1, "1M", false)[1]));
	}
	@Test
	public void isDiscordantPairMember() {
		assertFalse(SAMRecordUtil.isDiscordantPairMember(Read(0, 1, "1M")));
		assertFalse(SAMRecordUtil.isDiscordantPairMember(RP(1, 1, 1)[0]));
		assertFalse(SAMRecordUtil.isDiscordantPairMember(RP(1, 1, 1)[1]));
		assertFalse(SAMRecordUtil.isDiscordantPairMember(OEA(1, 1, "1M", true)[0]));
		assertFalse(SAMRecordUtil.isDiscordantPairMember(OEA(1, 1, "1M", false)[0]));
		assertFalse(SAMRecordUtil.isDiscordantPairMember(OEA(1, 1, "1M", true)[1]));
		assertFalse(SAMRecordUtil.isDiscordantPairMember(OEA(1, 1, "1M", false)[1]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", true, 0, 1, "1M", true)[0]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", true, 0, 1, "1M", false)[0]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", false, 0, 1, "1M", true)[0]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", false, 0, 1, "1M", false)[0]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", true, 0, 1, "1M", true)[1]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", true, 0, 1, "1M", false)[1]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", false, 0, 1, "1M", true)[1]));
		assertTrue(SAMRecordUtil.isDiscordantPairMember(DP(1, 1, "1M", false, 0, 1, "1M", false)[1]));
	}
	@Test
	public void isAnchoredPairMember() {
		assertFalse(SAMRecordUtil.isAnchoredPairMember(Read(0, 1, "1M")));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(RP(1, 1, 1)[0]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(RP(1, 1, 1)[1]));
		assertTrue(SAMRecordUtil.isAnchoredPairMember(OEA(1, 1, "1M", true)[0]));
		assertTrue(SAMRecordUtil.isAnchoredPairMember(OEA(1, 1, "1M", false)[0]));
		assertTrue(SAMRecordUtil.isAnchoredPairMember(OEA(1, 1, "1M", true)[1]));
		assertTrue(SAMRecordUtil.isAnchoredPairMember(OEA(1, 1, "1M", false)[1]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", true, 0, 1, "1M", true)[0]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", true, 0, 1, "1M", false)[0]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", false, 0, 1, "1M", true)[0]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", false, 0, 1, "1M", false)[0]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", true, 0, 1, "1M", true)[1]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", true, 0, 1, "1M", false)[1]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", false, 0, 1, "1M", true)[1]));
		assertFalse(SAMRecordUtil.isAnchoredPairMember(DP(1, 1, "1M", false, 0, 1, "1M", false)[1]));
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
}
