package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import net.sf.samtools.SAMRecord;

import org.junit.Test;

public class SoftClipEvidenceTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_b() {
		SAMRecord r = Read(0, 10, "10M5S");
		new SoftClipEvidence(BreakpointDirection.Backward, r);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_f() {
		SAMRecord r = Read(0, 10, "5S10M");
		new SoftClipEvidence(BreakpointDirection.Forward, r);
	}
	@Test
	public void GetBreakpointLocation_should_get_f_location() {
		SAMRecord r = Read(0, 10, "10M5S");
		BreakpointLocation l = new SoftClipEvidence(BreakpointDirection.Forward, r).getBreakpointLocation();
		assertEquals(BreakpointDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
	}
	@Test
	public void GetBreakpointLocation_should_get_interval_with_realigned_ff() {
		// MMMMSSS->
		//            SSS->  -ve strand flips from expected direction
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(true);
		BreakpointInterval l = (BreakpointInterval)new SoftClipEvidence(BreakpointDirection.Forward, r, realigned).getBreakpointLocation();
		assertEquals(BreakpointDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakpointDirection.Forward, l.direction2);
		assertEquals(12, l.start2);
		assertEquals(12, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakpointLocation_should_get_interval_with_realigned_fb() {
		// MMMMSSS->
		//            <-SSS
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(false);
		BreakpointInterval l = (BreakpointInterval)new SoftClipEvidence(BreakpointDirection.Forward, r, realigned).getBreakpointLocation();
		assertEquals(BreakpointDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakpointDirection.Backward, l.direction2);
		assertEquals(10, l.start2);
		assertEquals(10, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakpointLocation_should_get_interval_with_realigned_bf() {
		// SC going back so if we match our realigned read on the positive strand
		// then the join is at the end of the realigned read
		//      <-SSSMMMM
		// SSS->
		SAMRecord r = Read(0, 10, "6S10M");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(false);
		BreakpointInterval l = (BreakpointInterval)new SoftClipEvidence(BreakpointDirection.Backward, r, realigned).getBreakpointLocation();
		assertEquals(BreakpointDirection.Backward, l.direction);
		assertEquals(10, l.start);
		assertEquals(10, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakpointDirection.Forward, l.direction2);
		assertEquals(12, l.start2);
		assertEquals(12, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakpointLocation_should_get_interval_with_realigned_bb() {
		// SC going back so if we match our realigned read on the negative strand
		// then the join is at the start of the realigned read
		// since it was
		//      <-SSSMMMM
		// <-SSS  -ve strand flips from expected direction
		SAMRecord r = Read(0, 10, "6S10M");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(true);
		BreakpointInterval l = (BreakpointInterval)new SoftClipEvidence(BreakpointDirection.Backward, r, realigned).getBreakpointLocation();
		assertEquals(BreakpointDirection.Backward, l.direction);
		assertEquals(10, l.start);
		assertEquals(10, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakpointDirection.Backward, l.direction2);
		assertEquals(10, l.start2);
		assertEquals(10, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakpointLocation_should_get_location_with_unmapped_realigned() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Unmapped(15);
		assertFalse(new SoftClipEvidence(BreakpointDirection.Forward, r, realigned).getBreakpointLocation() instanceof BreakpointInterval);
	}
	@Test
	public void realigned_constructors_should_be_equivalent() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 15, "10M");
		SoftClipEvidence c1 = new SoftClipEvidence(BreakpointDirection.Forward, r, realigned);
		SoftClipEvidence c2 = new SoftClipEvidence(new SoftClipEvidence(BreakpointDirection.Forward, r), realigned);
		assertEquals(c1.getSAMRecord(), c2.getSAMRecord());
		assertEquals(c1.getSoftClipRealignmentSAMRecord(), c2.getSoftClipRealignmentSAMRecord());
		assertEquals(c1.getEvidenceID(), c2.getEvidenceID());
		assertEquals(c1.getBreakpointLocation().direction, c2.getBreakpointLocation().direction);
		assertEquals(c1.getBreakpointLocation().start, c2.getBreakpointLocation().start);
		assertEquals(c1.getBreakpointLocation().end, c2.getBreakpointLocation().end);
		assertEquals(c1.getBreakpointLocation().referenceIndex, c2.getBreakpointLocation().referenceIndex);
		assertEquals(((BreakpointInterval)c1.getBreakpointLocation()).direction2, ((BreakpointInterval)c2.getBreakpointLocation()).direction2);
		assertEquals(((BreakpointInterval)c1.getBreakpointLocation()).start2, ((BreakpointInterval)c2.getBreakpointLocation()).start2);
		assertEquals(((BreakpointInterval)c1.getBreakpointLocation()).end2, ((BreakpointInterval)c2.getBreakpointLocation()).end2);
		assertEquals(((BreakpointInterval)c1.getBreakpointLocation()).referenceIndex2, ((BreakpointInterval)c2.getBreakpointLocation()).referenceIndex2);
		assertEquals(c1.getBreakpointLocation().qual, c2.getBreakpointLocation().qual, 0);
	}
	@Test
	public void getEvidenceID_paired_should_encode_read_direction_and_pair_info() {
		SAMRecord r = RP(1, 1, 100)[0];
		r.setReadName("ReadName");
		r.setCigarString("1S2M3S");
		r.setFirstOfPairFlag(true);
		assertEquals("fReadName/1", new SoftClipEvidence(BreakpointDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName/1", new SoftClipEvidence(BreakpointDirection.Backward, r).getEvidenceID());
		r.setFirstOfPairFlag(false);
		assertEquals("fReadName/2", new SoftClipEvidence(BreakpointDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName/2", new SoftClipEvidence(BreakpointDirection.Backward, r).getEvidenceID());
	}
	@Test
	public void getEvidenceID_unpaired_should_encode_read_direction() {
		SAMRecord r = Read(1, 1, 100);
		r.setReadName("ReadName");
		r.setCigarString("1S2M3S");
		r.setReadName("ReadName");
		assertEquals("fReadName", new SoftClipEvidence(BreakpointDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName", new SoftClipEvidence(BreakpointDirection.Backward, r).getEvidenceID());
	}
	@Test
	public void unrealigned_should_have_decent_quality_metric() {
		assertFalse("TODO: get a decent quality metric", false);
	}
	@Test
	public void realigned_should_have_decent_quality_metric() {
		assertFalse("TODO: get a decent quality metric", false);
	}
	@Test
	public void realigned_unmapped_should_have_decent_quality_metric() {
		assertFalse("TODO: get a decent quality metric", false);
	}
}