package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

public class SAMRecordSummaryTest extends TestHelper {

	@Before
	public void setUp() throws Exception {
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_hard_clips() {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setCigarString("1H1S1M");
		assertTrue(SAMRecordUtil.isAlignmentSoftClipped(r));
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_empty_cigar() {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setCigarString("");
		assertFalse(SAMRecordUtil.isAlignmentSoftClipped(r));
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_null_cigar() {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setCigar(null);
		assertFalse(SAMRecordUtil.isAlignmentSoftClipped(r));
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_unmapped_reads() {
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setCigar(null);
		assertFalse(SAMRecordUtil.isAlignmentSoftClipped(r));
	}
	@Test
	public void getEndSoftClipLength_should_handle_hard_clips() {
		assertEquals(5, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "1H3S10M5S1H")));
	}
	@Test
	public void getStartSoftClipLength_should_handle_hard_clips() {
		assertEquals(3, SAMRecordUtil.getStartSoftClipLength(Read(0, 1, "1H3S10M5S1H")));
	}
	@Test
	public void getEndSoftClipLength_should_handle_soft_clips() {
		assertEquals(5, SAMRecordUtil.getEndSoftClipLength(Read(0, 1, "3S10M5S")));
	}
	@Test
	public void getStartSoftClipLength_should_handle_soft_clips() {
		assertEquals(3, SAMRecordUtil.getStartSoftClipLength(Read(0, 1, "3S10M5S")));
	}
}
