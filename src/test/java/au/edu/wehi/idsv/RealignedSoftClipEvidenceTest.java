package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;


public class RealignedSoftClipEvidenceTest extends TestHelper {
	@Test
	public void getUntemplatedSequenceLength_should_match_realign_soft_clip() {
		assertEquals(1, SoftClipEvidence.create(SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward, Read(0, 1, "4M6S")), Read(0, 1, "1S3M2S")).getUntemplatedSequence().length());
		assertEquals(2, SoftClipEvidence.create(SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward, Read(0, 1, "4M6S")), onNegative(Read(0, 1, "1S3M2S"))[0]).getUntemplatedSequence().length());
		assertEquals(2, SoftClipEvidence.create(SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Backward, Read(0, 1, "6S4M")), Read(0, 1, "1S3M2S")).getUntemplatedSequence().length());
		assertEquals(1, SoftClipEvidence.create(SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Backward, Read(0, 1, "6S4M")), onNegative(Read(0, 1, "1S3M2S"))[0]).getUntemplatedSequence().length());
	}
	@Test
	public void GetBreakendSummary_should_get_interval_with_realigned_ff() {
		// MMMMSSS->
		//            SSS->  -ve strand flips from expected direction
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(true);
		BreakpointSummary l = (BreakpointSummary)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward, r, realigned).getBreakendSummary();
		assertEquals(BreakendDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakendDirection.Forward, l.direction2);
		assertEquals(12, l.start2);
		assertEquals(12, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakendSummary_should_get_interval_with_realigned_fb() {
		// MMMMSSS->
		//            <-SSS
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(false);
		BreakpointSummary l = (BreakpointSummary)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward, r, realigned).getBreakendSummary();
		assertEquals(BreakendDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakendDirection.Backward, l.direction2);
		assertEquals(10, l.start2);
		assertEquals(10, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakendSummary_should_get_interval_with_realigned_bf() {
		// SC going back so if we match our realigned read on the positive strand
		// then the join is at the end of the realigned read
		//      <-SSSMMMM
		// SSS->
		SAMRecord r = Read(0, 10, "6S10M");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(false);
		BreakpointSummary l = (BreakpointSummary)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Backward, r, realigned).getBreakendSummary();
		assertEquals(BreakendDirection.Backward, l.direction);
		assertEquals(10, l.start);
		assertEquals(10, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakendDirection.Forward, l.direction2);
		assertEquals(12, l.start2);
		assertEquals(12, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void GetBreakendSummary_should_get_interval_with_realigned_bb() {
		// SC going back so if we match our realigned read on the negative strand
		// then the join is at the start of the realigned read
		// since it was
		//      <-SSSMMMM
		// <-SSS  -ve strand flips from expected direction
		SAMRecord r = Read(0, 10, "6S10M");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(true);
		BreakpointSummary l = (BreakpointSummary)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Backward, r, realigned).getBreakendSummary();
		assertEquals(BreakendDirection.Backward, l.direction);
		assertEquals(10, l.start);
		assertEquals(10, l.end);
		assertEquals(0, l.referenceIndex);
		assertEquals(BreakendDirection.Backward, l.direction2);
		assertEquals(10, l.start2);
		assertEquals(10, l.end2);
		assertEquals(1, l.referenceIndex2);
	}
	@Test
	public void should_set_realignment_attributes() {
		// MMMMSSS->
		//            SSS->  -ve strand flips from expected direction
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(true);
		SoftClipEvidence e = SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward, r, realigned);
		assertEquals(realigned.getReferenceIndex(), e.getSAMRecord().getIntegerAttribute("rr"));
		assertEquals(realigned.getAlignmentStart(), (int)e.getSAMRecord().getIntegerAttribute("rp"));
	}
	@Test
	public void getRemoteMapq_should_be_realigned_mapq() {
		assertEquals(5, ((RealignedSoftClipEvidence)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward,
				Read(0, 1, "4M6S"),
				withMapq(5, Read(0, 1, "1S3M2S"))[0]))
			.getRemoteMapq());
	}
	@Test
	public void getRemoteBaseLength_should_be_sc_length() {
		assertEquals(6, ((RealignedSoftClipEvidence)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward,
				Read(0, 1, "4M6S"),
				withMapq(5, Read(0, 1, "1S3M2S"))[0]))
			.getRemoteBaseLength());
	}
	@Test
	public void getRemoteBaseCount_should_be_sc_length() {
		assertEquals(6, ((RealignedSoftClipEvidence)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward,
				Read(0, 1, "4M6S"),
				withMapq(5, Read(0, 1, "1S3M2S"))[0]))
			.getRemoteBaseCount());
	}
	@Test
	public void getRemoteMaxBaseQual_realigned_mapped_quals() {
		assertEquals(4, ((RealignedSoftClipEvidence)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward,
				Read(0, 1, "4M6S"),
				withQual(new byte[] { 1,2,3,4,5,6}, Read(0, 1, "1S3M2S"))[0]))
			.getRemoteMaxBaseQual());
	}
	@Test
	public void getRemoteTotalBaseQual_realigned_mapped_quals() {
		assertEquals(2+3+4, ((RealignedSoftClipEvidence)SoftClipEvidence.create(getContext(), SES(), BreakendDirection.Forward,
				Read(0, 1, "4M6S"),
				withQual(new byte[] { 1,2,3,4,5,6}, Read(0, 1, "1S3M2S"))[0]))
			.getRemoteTotalBaseQual());
	}

}
