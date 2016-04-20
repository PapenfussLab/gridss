package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.io.File;

import org.apache.commons.configuration.ConfigurationException;
import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;

public class SoftClipEvidenceTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_b() {
		SAMRecord r = Read(0, 10, "10M5S");
		SoftClipEvidence.create(SES(), BreakendDirection.Backward, r);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_f() {
		SAMRecord r = Read(0, 10, "5S10M");
		SoftClipEvidence.create(SES(), BreakendDirection.Forward, r);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_cigar() {
		SoftClipEvidence.create(SES(), BreakendDirection.Backward, new SAMRecord(getContext().getBasicSamHeader()));
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_read_bases() {
		SoftClipEvidence.create(SES(), BreakendDirection.Backward, withSequence((byte[])null, Read(0, 1, "1S1M"))[0]);
	}
	@Test
	public void GetBreakendSummary_should_get_f_location() {
		SAMRecord r = Read(0, 10, "10M5S");
		BreakendSummary l = SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getBreakendSummary();
		assertEquals(BreakendDirection.Forward, l.direction);
		assertEquals(19, l.start);
		assertEquals(19, l.end);
		assertEquals(0, l.referenceIndex);
	}
	@Test
	public void GetBreakendSummary_should_get_b_location() {
		SAMRecord r = Read(0, 10, "5S10M");
		BreakendSummary l = SoftClipEvidence.create(SES(), BWD, r).getBreakendSummary();
		assertEquals(BWD, l.direction);
		assertEquals(10, l.start);
		assertEquals(10, l.end);
		assertEquals(0, l.referenceIndex);
	}
	@Test
	public void GetBreakendSummary_should_get_location_with_unmapped_realigned() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Unmapped(15);
		assertFalse(SoftClipEvidence.create(SES(), BreakendDirection.Forward, r, realigned).getBreakendSummary() instanceof BreakpointSummary);
	}
	@Test
	public void realigned_constructors_should_be_equivalent() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 15, "10M");
		SoftClipEvidence c1 = SoftClipEvidence.create(SES(), BreakendDirection.Forward, r, realigned);
		SoftClipEvidence c2 = SoftClipEvidence.create(SoftClipEvidence.create(SES(), BreakendDirection.Forward, r), realigned);
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
		assertEquals("fReadName\\1", SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName\\1", SoftClipEvidence.create(SES(), BreakendDirection.Backward, r).getEvidenceID());
		r.setFirstOfPairFlag(false);
		assertEquals("fReadName\\2", SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName\\2", SoftClipEvidence.create(SES(), BreakendDirection.Backward, r).getEvidenceID());
	}
	@Test
	public void getEvidenceID_unpaired_should_encode_read_direction() {
		SAMRecord r = Read(1, 1, 100);
		r.setReadName("ReadName");
		r.setCigarString("1S2M3S");
		r.setReadName("ReadName");
		assertEquals("fReadName", SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID());
		assertEquals("bReadName", SoftClipEvidence.create(SES(), BreakendDirection.Backward, r).getEvidenceID());
	}
	@Test
	public void realigned_consider_multimapping_realignment_unmapped() {
		SAMRecord r = Read(0, 10, "10M5S");
		SAMRecord realigned = Read(1, 15, "10M");
		realigned.setMappingQuality(1);
		SoftClipEvidence e = SoftClipEvidence.create(SES(), BreakendDirection.Forward, r, realigned);
		assertFalse(e instanceof RealignedSoftClipEvidence);
		assertFalse(e.getBreakendSummary() instanceof BreakpointSummary);
	}
	@Ignore
	@Test
	public void unrealigned_should_have_decent_quality_metric() {
		assertTrue("TODO: get a decent quality metric", false);
	}
	@Ignore
	@Test
	public void realigned_should_have_decent_quality_metric() {
		assertTrue("TODO: get a decent quality metric", false);
	}
	@Ignore
	@Test
	public void realigned_unmapped_should_have_decent_quality_metric() {
		assertTrue("TODO: get a decent quality metric", false);
	}
	@Test
	public void getAverageClipQuality_should_average_clip_quality_bases() {
		SAMRecord r = Read(0, 1, "2M3S");
		r.setBaseQualities(new byte[] { 40, 30, 20, 10, 0 });
		SoftClipEvidence e = SoftClipEvidence.create(SES(), BreakendDirection.Forward, r);
		assertEquals(10, e.getAverageClipQuality(), 0);
	}
	@Test
	public void getAverageClipQuality_should_return_0_if_unknown() {
		SAMRecord r = Read(0, 1, "2M3S");
		r.setBaseQualities(null);
		assertEquals(0, SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getAverageClipQuality(), 0);
		
		r = Read(0, 1, "2M3S");
		r.setBaseQualities(SAMRecord.NULL_QUALS);
		assertEquals(0, SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getAverageClipQuality(), 0);
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip() {
		SoftClipEvidence.create(SES(), BreakendDirection.Forward, new SAMRecord(getHeader()));
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_forward() {
		SoftClipEvidence.create(SES(), BreakendDirection.Forward, Read(0, 1, "1S10M")); 
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_backward() {
		SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "10M1S")); 
	}
	@Test
	public void getSoftClipLength_should_return_clip_length() {
		assertEquals(2, SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "2S5M1S")).getSoftClipLength());
		assertEquals(2, SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "1H2S5M1S10H")).getSoftClipLength());
	}
	@Test
	public void getUntemplatedSequenceLength_should_match_soft_clip_with_no_realign() {
		SoftClipEvidence sc = SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "1S4M"));
		assertEquals(sc.getSoftClipLength(), sc.getBreakendSequence().length);
		sc = SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "1S4M2S"));
		assertEquals(sc.getSoftClipLength(), sc.getBreakendSequence().length);
		sc = SoftClipEvidence.create(SES(), BreakendDirection.Forward, Read(0, 1, "4M1S"));
		assertEquals(sc.getSoftClipLength(), sc.getBreakendSequence().length);
		sc = SoftClipEvidence.create(SES(), BreakendDirection.Forward, Read(0, 1, "1S4M5S"));
		assertEquals(sc.getSoftClipLength(), sc.getBreakendSequence().length);
	}
	@Test
	public void getUntemplatedSequenceLength_should_match_soft_clip_with_realign_unmapped() {
		assertEquals(1, SoftClipEvidence.create(SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "1S4M")), Unmapped(1)).getBreakendSequence().length);
	}
	@Test
	public void getLocalMapq_should_be_mapq() {
		assertEquals(5, SCE(FWD, withMapq(5, Read(0, 1, "4M2S"))).getLocalMapq());
	}
	@Test
	public void getLocalBaseLength_should_exclude_sc_bases() {
		assertEquals(4, SCE(FWD, withMapq(5, Read(0, 1, "4M2S"))).getLocalBaseLength());
	}
	@Test
	public void getLocalMaxBaseQual_should_be_reference_quals() {
		assertEquals(4, SCE(FWD, withQual(new byte[] {1,2,3,4,5,6}, Read(0, 1, "4M2S"))).getLocalMaxBaseQual());
	}
	@Test
	public void getLocalTotalBaseQual_should_be_reference_quals() {
		assertEquals(1+2+3+4, SCE(FWD, withQual(new byte[] {1,2,3,4,5,6}, Read(0, 1, "4M2S"))).getLocalTotalBaseQual());
	}
	private MockSAMEvidenceSource permissiveSES() {
		MockSAMEvidenceSource ses = SES();
		ses.getContext().getConfig().getSoftClip().minAnchorIdentity = 0;
		ses.getContext().getConfig().getSoftClip().minLength = 1;
		ses.getContext().getConfig().minMapq = 0;
		ses.getContext().getConfig().minAnchorShannonEntropy = 0;
		ses.getContext().getConfig().adapters = new AdapterHelper(new String[0]);
		return ses;
	}
	private void adapter(String seq, int scLen, boolean keep) {
		SAMEvidenceSource ses = permissiveSES();
		try {
			ses.getContext().getConfig().adapters = new GridssConfiguration().adapters;
		} catch (ConfigurationException e) {
		}
		SAMRecord r = Read(0, 1000, String.format("%dM%dS", seq.length() - scLen, scLen));
		r.setReadBases(B(seq));
		r.setReadNegativeStrandFlag(false);
		assertEquals(keep, new SoftClipEvidence(ses, FWD, r).meetsEvidenceCritera());
		// reverse comp
		r = Read(0, 1000, String.format("%dS%dM", scLen, seq.length() - scLen));
		r.setReadBases(B(SequenceUtil.reverseComplement(seq)));
		r.setReadNegativeStrandFlag(true);
		assertEquals(keep, new SoftClipEvidence(ses, BWD, r).meetsEvidenceCritera());
	}
	@Test
	public void should_filter_adapter_sequence_IlluminaUniversalAdapter_AGATCGGAAGAG() {
		adapter("TTTTTAGATCGGAAGAG", 12, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 24, false);
		adapter("AAACATTTTGCCATTTTTATGGGCAAAAATGATAATTTCTTGTTAATTTAATTTGTATTTTTCCCTGAATGAGGTAAGCTTTTAAATACTTATATTAGATC", 5, false);
		// adapter sequence on the wrong end of the read
		adapter("TCTTCCGATCTTGTTGTTTTCTTAGGGAAAAAATTCTAGCAGTGGTTTTGCTTGAAATAAGAATATGGATCCTGTAGGCTTTTGATACACCTTGCTAAATT", 10, true);
	}
	@Test
	public void should_filter_adapter_sequence_IlluminaSmallRNAAdapter_ATGGAATTCTCG() {
		adapter("TTTTTATGGAATTCTCG", 12, false);
	}
	@Test
	public void should_filter_adapter_sequence_NexteraTransposaseSequence_CTGTCTCTTATA() {
		adapter("TTTTTCTGTCTCTTATA", 12, false);
	}
	@Test
	public void should_filter_adapter_short_read_short_sc() {
		// only 1 BP mapped
		adapter("TAGATCG", 6, false);
	}
	@Test
	public void should_filter_adapter_short_anchor_long_clip() {
		// only 1 BP mapped
		adapter("TAGATCGGAAGAGTATCATCTACTACTATCTACTATCATCTACTATCTA", 48, false);
		adapter("TTGATCGGAAGAGTATCATCTACTACTATCTACTATCATCTACTATCTA", 48, true);
	}
	@Test
	public void should_filter_adapter_long_anchor_long_sc() {
		// only 1 BP mapped
		adapter("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAGATCGGAAGAGTATCATCTACTACTATCTACTATCATCTACTATCTA", 48, false);
		adapter("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTGATCGGAAGAGTATCATCTACTACTATCTACTATCATCTACTATCTA", 48, true);
	}
	@Test
	public void should_allow_adapter_reference_microhomology_up_to_6bp() {
		// adapter bases could be mapped - we need to check upstream
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 23, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 22, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 21, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 20, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 19, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 18, false);
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCGGAAGAGCGTCGTGTAGGG", 17, true);
	}
	@Test
	public void n_should_match_adapter_base() {
		adapter("GTTAATTTAATTTGTATTTTTCCCTGAATTAGGTAGGCTTTTAAATACTTATATTGCAATCTGTATTTCATGTTTGTAGATCNGAAGAGCGTCGTGTAGGG", 24, false);
		//                                                                                         ^
	}
	@Test
	public void should_filter_dovetailing_reads() {
		SAMEvidenceSource ses = permissiveSES();
		SAMRecord[] rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10M10S");
		rp[1].setCigarString("10S10M");
		assertFalse(new SoftClipEvidence(ses, FWD, rp[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, BWD, rp[1]).meetsEvidenceCritera());
		
		// looks like a dovetail, but is not
		rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10S10M"); // <-- definitely keep this one
		rp[1].setCigarString("10S10M"); // ideally keep this one too, but we need to know about the mate cigar for that
		assertTrue(new SoftClipEvidence(SES(), BWD, rp[0]).meetsEvidenceCritera());
		//assertTrue(new SoftClipEvidence(getContext(), ses, BWD, rp[1]).meetsEvidenceCritera(scp));
	}
	@Test
	public void should_allow_4_bp_dovetail_margin() {
		SAMEvidenceSource ses = permissiveSES();
		for (int i = 95 ; i <= 105; i++) { // -3bp -> +3bp, only +-3bp pass dovetail filter  
			SAMRecord[] rp = RP(0, 100, i, 20);
			rp[0].setCigarString("10M10S");
			rp[1].setCigarString("10S10M");
			assertEquals(i == 95 || i == 105, new SoftClipEvidence(ses, FWD, rp[0]).meetsEvidenceCritera());
			assertEquals(i == 95 || i == 105, new SoftClipEvidence(ses, BWD, rp[1]).meetsEvidenceCritera());
		}
	}
	@Test
	public void should_filter_dovetailing_reads_with_metrics() {
		MockSAMEvidenceSource ses = permissiveSES();
		ses.metrics = new IdsvSamFileMetrics(
				new File("src/test/resources/testmetrics.metrics.idsv.txt"),
				new File("src/test/resources/testmetrics.metrics.insertsize.txt"),
				new File("src/test/resources/testmetrics.metrics.mapq.txt"),
				new File("src/test/resources/testmetrics.metrics.cigar.txt"));
		SAMRecord[] rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10M10S");
		rp[1].setCigarString("10S10M");
		assertFalse(new SoftClipEvidence(ses, FWD, rp[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, BWD, rp[1]).meetsEvidenceCritera());
		
		// looks like a dovetail, but is not
		rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10S10M"); // <-- definitely keep this one
		rp[1].setCigarString("10S10M"); // ideally keep this one too, but we need to know about the mate cigar for that
		assertTrue(new SoftClipEvidence(ses, BWD, rp[0]).meetsEvidenceCritera());
		//assertTrue(new SoftClipEvidence(getContext(), ses, BWD, rp[1]).meetsEvidenceCritera(scp));
	}
	@Test
	public void isNotExact() {
		assertTrue(SCE(FWD, Read(0, 1, "1M1S")).isBreakendExact());
	}
	@Test
	public void getAnchorSequence_should_return_all_bases_not_in_soft_clip() {
		assertEquals("ACGT", S(SCE(FWD, withSequence("ACGTNN", Read(0, 1, "1S3M2S"))[0]).getAnchorSequence()));
		assertEquals("CGTNN", S(SCE(BWD, withSequence("ACGTNN", Read(0, 1, "1S3M2S"))[0]).getAnchorSequence()));
	}
	@Test
	public void should_filter_low_complexity_anchors() {
		SAMEvidenceSource ses = permissiveSES();
		ses.getContext().getConfig().minAnchorShannonEntropy = 0.5;
		
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AT", Read(0, 1, "1M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAT", Read(0, 1, "2M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAT", Read(0, 1, "3M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAT", Read(0, 1, "4M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAT", Read(0, 1, "5M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAT", Read(0, 1, "6M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAAT", Read(0, 1, "7M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAAAT", Read(0, 1, "8M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAAAAAT", Read(0, 1, "10M1S"))[0]).meetsEvidenceCritera());
		
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("GT", Read(0, 1, "1M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAT", Read(0, 1, "2M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAT", Read(0, 1, "3M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAAT", Read(0, 1, "4M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAAAT", Read(0, 1, "5M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAAAAT", Read(0, 1, "6M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAAAAAT", Read(0, 1, "7M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAAAAAAT", Read(0, 1, "8M1S"))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("GAAAAAAAAT", Read(0, 1, "9M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("GAAAAAAAAAT", Read(0, 1, "10M1S"))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("GAAAAAAAAAAT", Read(0, 1, "11M1S"))[0]).meetsEvidenceCritera());
		
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "40M1S"))[0]).meetsEvidenceCritera()); // 37	1	1	1	0.503183732
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "41M1S"))[0]).meetsEvidenceCritera()); // 38	1	1	1	0.493619187
		
		assertTrue(new SoftClipEvidence(ses, FWD, withSequence("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "1S40M1S"))[0]).meetsEvidenceCritera()); // 37	1	1	1	0.503183732
		assertFalse(new SoftClipEvidence(ses, FWD, withSequence("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "1S41M1S"))[0]).meetsEvidenceCritera()); // 38	1	1	1	0.493619187
	}
	@Test
	public void should_filter_low_base_quality() {
		SAMEvidenceSource ses = permissiveSES();
		ses.getContext().getConfig().getSoftClip().minAverageQual = 5;
		
		assertTrue(new SoftClipEvidence(ses, FWD, withQual(new byte[] { 6, 6, }, withSequence("AT", Read(0, 1, "1M1S")))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withQual(new byte[] { 6, 5, }, withSequence("AT", Read(0, 1, "1M1S")))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withQual(new byte[] { 6, 4, }, withSequence("AT", Read(0, 1, "1M1S")))[0]).meetsEvidenceCritera());
		assertTrue(new SoftClipEvidence(ses, FWD, withQual(new byte[] { 1, 4, 6, }, withSequence("ATT", Read(0, 1, "1M2S")))[0]).meetsEvidenceCritera());
		assertFalse(new SoftClipEvidence(ses, FWD, withQual(new byte[] { 1, 4, 5, }, withSequence("ATT", Read(0, 1, "1M2S")))[0]).meetsEvidenceCritera());
	}
}
