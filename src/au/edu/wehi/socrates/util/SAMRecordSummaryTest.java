package au.edu.wehi.socrates.util;

import static org.junit.Assert.*;

import net.sf.samtools.SAMRecord;

import org.junit.Before;
import org.junit.Test;

public class SAMRecordSummaryTest {

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
}
