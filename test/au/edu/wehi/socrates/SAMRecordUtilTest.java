package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.SAMRecord;

import org.junit.Test;


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
		SAMRecordUtil.ensureNmTag(null, r);
	}

	@Test
	public void ensureNmTag_should_set_NM() {
		// AAAA
		// ACGT
		SAMRecord r = Read(1, 1, "4M");
		SAMRecordUtil.ensureNmTag(new ReferenceSequenceFileWalker(SMALL_FA), r);
		assertEquals(3, (int)r.getIntegerAttribute("NM"));
	}
	public void ensureNmTag_should_ignore_unmapped() {
		SAMRecord r = Unmapped(5);
		SAMRecordUtil.ensureNmTag(new ReferenceSequenceFileWalker(SMALL_FA), r);
		assertNull(r.getIntegerAttribute("NM"));
	}
}
