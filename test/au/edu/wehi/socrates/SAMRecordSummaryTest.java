package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import net.sf.samtools.SAMRecord;

import org.junit.Before;
import org.junit.Test;

import au.edu.wehi.socrates.util.SAMRecordSummary;

public class SAMRecordSummaryTest extends TestHelper {

	@Before
	public void setUp() throws Exception {
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_hard_clips() {
		SAMRecord r = new SAMRecord(null);
		r.setCigarString("1H1S1M");
		assertTrue(SAMRecordSummary.isAlignmentSoftClipped(r));
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_empty_cigar() {
		SAMRecord r = new SAMRecord(null);
		r.setCigarString("");
		assertFalse(SAMRecordSummary.isAlignmentSoftClipped(r));
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_null_cigar() {
		SAMRecord r = new SAMRecord(null);
		r.setCigar(null);
		assertFalse(SAMRecordSummary.isAlignmentSoftClipped(r));
	}
	@Test
	public void isAlignmentSoftClipped_should_handle_unmapped_reads() {
		SAMRecord r = new SAMRecord(null);
		r.setCigar(null);
		assertFalse(SAMRecordSummary.isAlignmentSoftClipped(r));
	}
	@Test
	public void getEndSoftClipLength_should_handle_hard_clips() {
		assertEquals(5, SAMRecordSummary.getEndSoftClipLength(Read(0, 1, "1H3S10M5S1H")));
	}
	@Test
	public void getStartSoftClipLength_should_handle_hard_clips() {
		assertEquals(3, SAMRecordSummary.getStartSoftClipLength(Read(0, 1, "1H3S10M5S1H")));
	}
	@Test
	public void getEndSoftClipLength_should_handle_soft_clips() {
		assertEquals(5, SAMRecordSummary.getEndSoftClipLength(Read(0, 1, "3S10M5S")));
	}
	@Test
	public void getStartSoftClipLength_should_handle_soft_clips() {
		assertEquals(3, SAMRecordSummary.getStartSoftClipLength(Read(0, 1, "3S10M5S")));
	}
}
