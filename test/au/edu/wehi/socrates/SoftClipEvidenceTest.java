package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import au.edu.wehi.socrates.vcf.EvidenceAttributes;

public class SoftClipEvidenceTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_b() {
		SAMRecord r = Read(0, 10, "10M5S");
		new SoftClipEvidence(getContext(), BreakendDirection.Backward, r);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_f() {
		SAMRecord r = Read(0, 10, "5S10M");
		new SoftClipEvidence(getContext(), BreakendDirection.Forward, r);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_cigar() {
		new SoftClipEvidence(getContext(), BreakendDirection.Backward, new SAMRecord(null));
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_read_bases() {
		new SoftClipEvidence(getContext(), BreakendDirection.Backward, withSequence((byte[])null, Read(0, 1, "1S1M"))[0]);
	}
	@Test
	public void GetBreakendSummary_should_get_f_location() {
		SAMRecord r = Read(0, 10, "10M5S");
		BreakendSummary l = new SoftClipEvidence(getContext(), BreakendDirection.Forward, r).getBreakendSummary();
		assertEquals(BreakendDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
	}
	@Test
	public void GetBreakendSummary_should_get_interval_with_realigned_ff() {
		// MMMMSSS->
		//            SSS->  -ve strand flips from expected direction
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 10, "1S3M2S");
		realigned.setReadNegativeStrandFlag(true);
		BreakpointSummary l = (BreakpointSummary)new SoftClipEvidence(getContext(), BreakendDirection.Forward, r, realigned).getBreakendSummary();
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
		BreakpointSummary l = (BreakpointSummary)new SoftClipEvidence(getContext(), BreakendDirection.Forward, r, realigned).getBreakendSummary();
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
		BreakpointSummary l = (BreakpointSummary)new SoftClipEvidence(getContext(), BreakendDirection.Backward, r, realigned).getBreakendSummary();
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
		BreakpointSummary l = (BreakpointSummary)new SoftClipEvidence(getContext(), BreakendDirection.Backward, r, realigned).getBreakendSummary();
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
	public void GetBreakendSummary_should_get_location_with_unmapped_realigned() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Unmapped(15);
		assertFalse(new SoftClipEvidence(getContext(), BreakendDirection.Forward, r, realigned).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void realigned_constructors_should_be_equivalent() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 15, "10M");
		SoftClipEvidence c1 = new SoftClipEvidence(getContext(), BreakendDirection.Forward, r, realigned);
		SoftClipEvidence c2 = new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Forward, r), realigned);
		assertEquals(c1.getSAMRecord(), c2.getSAMRecord());
		assertEquals(c1.getEvidenceID(), c2.getEvidenceID());
		assertEquals(c1.getBreakendSummary().direction, c2.getBreakendSummary().direction);
		assertEquals(c1.getBreakendSummary().start, c2.getBreakendSummary().start);
		assertEquals(c1.getBreakendSummary().end, c2.getBreakendSummary().end);
		assertEquals(c1.getBreakendSummary().referenceIndex, c2.getBreakendSummary().referenceIndex);
		assertEquals(((BreakpointSummary)c1.getBreakendSummary()).direction2, ((BreakpointSummary)c2.getBreakendSummary()).direction2);
		assertEquals(((BreakpointSummary)c1.getBreakendSummary()).start2, ((BreakpointSummary)c2.getBreakendSummary()).start2);
		assertEquals(((BreakpointSummary)c1.getBreakendSummary()).end2, ((BreakpointSummary)c2.getBreakendSummary()).end2);
		assertEquals(((BreakpointSummary)c1.getBreakendSummary()).referenceIndex2, ((BreakpointSummary)c2.getBreakendSummary()).referenceIndex2);
	}
	@Test
	public void getEvidenceID_paired_should_encode_read_direction_and_pair_info() {
		SAMRecord r = RP(1, 1, 100)[0];
		r.setReadName("ReadName");
		r.setCigarString("1S2M3S");
		r.setFirstOfPairFlag(true);
		assertEquals("fReadName/1", new SoftClipEvidence(getContext(), BreakendDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName/1", new SoftClipEvidence(getContext(), BreakendDirection.Backward, r).getEvidenceID());
		r.setFirstOfPairFlag(false);
		assertEquals("fReadName/2", new SoftClipEvidence(getContext(), BreakendDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName/2", new SoftClipEvidence(getContext(), BreakendDirection.Backward, r).getEvidenceID());
	}
	@Test
	public void getEvidenceID_unpaired_should_encode_read_direction() {
		SAMRecord r = Read(1, 1, 100);
		r.setReadName("ReadName");
		r.setCigarString("1S2M3S");
		r.setReadName("ReadName");
		assertEquals("fReadName", new SoftClipEvidence(getContext(), BreakendDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName", new SoftClipEvidence(getContext(), BreakendDirection.Backward, r).getEvidenceID());
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
	@Test
	public void getAverageClipQuality_should_average_clip_quality_bases() {
		SAMRecord r = Read(0, 1, "2M3S");
		r.setBaseQualities(new byte[] { 40, 30, 20, 10, 0 });
		SoftClipEvidence e = new SoftClipEvidence(getContext(), BreakendDirection.Forward, r);
		assertEquals(10, e.getAverageClipQuality(), 0);
	}
	@Test
	public void getAverageClipQuality_should_return_0_if_unknown() {
		SAMRecord r = Read(0, 1, "2M3S");
		r.setBaseQualities(null);
		assertEquals(0, new SoftClipEvidence(getContext(), BreakendDirection.Forward, r).getAverageClipQuality(), 0);
		
		r = Read(0, 1, "2M3S");
		r.setBaseQualities(SAMRecord.NULL_QUALS);
		assertEquals(0, new SoftClipEvidence(getContext(), BreakendDirection.Forward, r).getAverageClipQuality(), 0);
	}
	@Test
	public void getAlignedPercentIdentity_should_match_only_mapped_bases() {
		assertEquals(100, new SoftClipEvidence(getContext(), BreakendDirection.Backward, withNM(withSequence("NAAAAA", Read(0, 1, "1S5M")))[0]).getAlignedPercentIdentity(), 0);
		assertEquals(100, new SoftClipEvidence(getContext(), BreakendDirection.Backward, withNM(withSequence("NTAAAT", Read(0, 1, "2S3M1S")))[0]).getAlignedPercentIdentity(), 0);
		assertEquals(100, new SoftClipEvidence(getContext(), BreakendDirection.Backward, withNM(withSequence("NATTTA", Read(0, 1, "1S1M3I1M")))[0]).getAlignedPercentIdentity(), 0);
		assertEquals(100, new SoftClipEvidence(getContext(), BreakendDirection.Backward, withNM(withSequence("NAGTAC", Read(1, 1, "1S1M1D4M")))[0]).getAlignedPercentIdentity(), 0);
		assertEquals(50, new SoftClipEvidence(getContext(), BreakendDirection.Backward, withNM(withSequence("NAATT", Read(0, 1, "1S4M")))[0]).getAlignedPercentIdentity(), 0);
		assertEquals(0, new SoftClipEvidence(getContext(), BreakendDirection.Backward, withNM(withSequence("ACCCCA", Read(0, 1, "1S4M1S")))[0]).getAlignedPercentIdentity(), 0);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip() {
		new SoftClipEvidence(getContext(), BreakendDirection.Forward, new SAMRecord(getHeader()));
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_forward() {
		new SoftClipEvidence(getContext(), BreakendDirection.Forward, Read(0, 1, "1S10M")); 
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_backward() {
		new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "10M1S")); 
	}
	@Test
	public void getSoftClipLength_should_return_clip_length() {
		assertEquals(2, new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "2S5M1S")).getSoftClipLength());
		assertEquals(2, new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "1H2S5M1S10H")).getSoftClipLength());
	}
	@Test
	public void getUntemplatedSequenceLength_should_match_soft_clip_with_no_realign() {
		SoftClipEvidence sc = new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "1S4M"));
		assertEquals(sc.getSoftClipLength(), sc.getUntemplatedSequenceLength());
		sc = new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "1S4M2S"));
		assertEquals(sc.getSoftClipLength(), sc.getUntemplatedSequenceLength());
		sc = new SoftClipEvidence(getContext(), BreakendDirection.Forward, Read(0, 1, "4M1S"));
		assertEquals(sc.getSoftClipLength(), sc.getUntemplatedSequenceLength());
		sc = new SoftClipEvidence(getContext(), BreakendDirection.Forward, Read(0, 1, "1S4M5S"));
		assertEquals(sc.getSoftClipLength(), sc.getUntemplatedSequenceLength());
	}
	@Test
	public void getUntemplatedSequenceLength_should_match_soft_clip_with_realign_unmapped() {
		assertEquals(1, new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "1S4M")), Unmapped(1)).getUntemplatedSequenceLength());
	}
	@Test
	public void getUntemplatedSequenceLength_should_match_realign_soft_clip() {
		assertEquals(1, new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Forward, Read(0, 1, "4M6S")), Read(0, 1, "1S3M2S")).getUntemplatedSequenceLength());
		assertEquals(2, new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Forward, Read(0, 1, "4M6S")), onNegative(Read(0, 1, "1S3M2S"))[0]).getUntemplatedSequenceLength());
		assertEquals(2, new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "6S4M")), Read(0, 1, "1S3M2S")).getUntemplatedSequenceLength());
		assertEquals(1, new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "6S4M")), onNegative(Read(0, 1, "1S3M2S"))[0]).getUntemplatedSequenceLength());
	}
	@Test
	public void should_set_sc_evidence() {
		EvidenceMetrics e = new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "3S4M")).getBreakendSummary().evidence;
		assertEquals(1, e.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT), 1);
		assertEquals(1, e.get(EvidenceAttributes.SOFT_CLIP_TOTAL_LENGTH), 3);
	}
}
