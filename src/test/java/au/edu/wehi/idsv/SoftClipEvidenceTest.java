package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.stream.Stream;

import org.apache.commons.configuration.ConfigurationException;
import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class SoftClipEvidenceTest extends TestHelper {
	static {
		//  1234567890
		//  ATGTGGC
		//  SSMMSSSH
		SAMRecord r = Read(2, 3, "2S2M3S5H");
		r.setReadBases(B("ATGTGGC"));
		r.setBaseQualities(B("1234567"));
		r.setMappingQuality(40);
		fExample = SoftClipEvidence.create(SES(), FWD, r);
		bExample = SoftClipEvidence.create(SES(), BWD, r);
	}
	private static final SoftClipEvidence fExample;
	private static final SoftClipEvidence bExample;
	@Test
	public void should_break_at_clip() {
		assertEquals(new BreakendSummary(2, FWD, 4), fExample.getBreakendSummary());
		assertEquals(new BreakendSummary(2, BWD, 3), bExample.getBreakendSummary());
	}
	@Test
	public void indels_should_have_unique_evidenceID() {
		assertNotEquals(fExample.getEvidenceID(), bExample.getEvidenceID());
	}
	@Test
	public void getBreakendSequence_should_return_sc() {
		assertEquals("GGC", S(fExample.getBreakendSequence()));
		assertEquals("AT", S(bExample.getBreakendSequence()));
	}
	
	@Test
	public void getBreakendQuality_should_return_sc_qual() {
		assertEquals("567", S(fExample.getBreakendQuality()));
		assertEquals("12", S(bExample.getBreakendQuality()));
	}
	
	@Test
	public void getAnchorSequence_should_exclude_soft_clipped_bases_on_other_side() {
		assertEquals("GT", S(fExample.getAnchorSequence()));
		assertEquals("GT", S(bExample.getAnchorSequence()));
	}
	
	@Test
	public void getAnchorQuality_should_return_local() {
		assertEquals("34", S(fExample.getAnchorQuality()));
		assertEquals("34", S(bExample.getAnchorQuality()));
	}
	
	@Test
	public void getMapq_should_match_read_mapq() {
		assertEquals(40, fExample.getLocalMapq());
		assertEquals(40, bExample.getLocalMapq());
	}
	
	@Test
	public void isBreakendExact() {
		assertTrue(fExample.isBreakendExact());
		assertTrue(bExample.isBreakendExact());
	}
	@Test
	public void getUntemplatedSequence_should_return_insert() {
		assertEquals(S(fExample.getBreakendSequence()), fExample.getUntemplatedSequence());
		assertEquals(S(bExample.getBreakendSequence()), bExample.getUntemplatedSequence());
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_f() {
		SoftClipEvidence.create(SES(), FWD, Read(0, 10, "5S10M"));
	}
	@Test(expected=IllegalArgumentException.class)
	public void constructor_should_require_soft_clip_b() {
		SoftClipEvidence.create(SES(), BWD, Read(0, 10, "10M5S"));
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
	public void getEvidenceID_paired_should_encode_pair_info() {
		SAMRecord r = RP(1, 1, 100)[0];
		r.setReadName("ReadName");
		r.setCigarString("10M10S");
		r.setReadPairedFlag(true);
		r.setFirstOfPairFlag(true);
		String r1 = SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID();
		r.setFirstOfPairFlag(false);
		r.setSecondOfPairFlag(true);
		String r2 = SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID();
		assertNotEquals(r1, r2);
	}
	@Test
	public void getEvidenceID_unpaired_should_encode_read_direction() {
		SAMRecord r = Read(1, 1, 100);
		r.setReadName("ReadName");
		r.setCigarString("5S5M5S");
		r.setReadName("ReadName");
		assertNotEquals(SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID(),
				SoftClipEvidence.create(SES(), BreakendDirection.Backward, r).getEvidenceID());
	}
	@Test
	public void should_recognise_1XS_as_unanchored_sequence() {
		SoftClipEvidence e = SoftClipEvidence.create(SES(), FWD, withSequence("NCGT", Read(0, 100, "1X3S"))[0]);
		assertFalse(e.isBreakendExact());
		assertEquals("CGT", S(e.getBreakendSequence()));
		assertEquals("", S(e.getAnchorSequence()));
		assertEquals("CGT", e.getUntemplatedSequence());
		assertEquals(new BreakendSummary(0, FWD, 100), e.getBreakendSummary());
	}
	@Test
	public void should_recognise_2XS_as_unanchored_sequence() {
		SoftClipEvidence e = SoftClipEvidence.create(SES(), FWD, withSequence("NNCGT", Read(0, 100, "2X3S"))[0]);
		assertFalse(e.isBreakendExact());
		assertEquals("CGT", S(e.getBreakendSequence()));
		assertEquals("", S(e.getAnchorSequence()));
		assertEquals(new BreakendSummary(0, FWD, 101, 100, 101), e.getBreakendSummary());
	}
	@Test
	public void should_recognise_XNX_as_unanchored_breakend_interval() {
		SoftClipEvidence e = SoftClipEvidence.create(SES(), FWD, withSequence("NNCGT", Read(0, 100, "1X1N1X3S"))[0]);
		assertFalse(e.isBreakendExact());
		assertEquals("CGT", S(e.getBreakendSequence()));
		assertEquals("", S(e.getAnchorSequence()));
		assertEquals(new BreakendSummary(0, FWD, 101, 100, 102), e.getBreakendSummary());
	}
	@Test
	public void should_recognise_SXNX_as_unanchored_sequence() {
		SoftClipEvidence e = SoftClipEvidence.create(SES(), BWD, withSequence("CGTNN", Read(0, 100, "3S1X2N1X"))[0]);
		assertFalse(e.isBreakendExact());
		assertEquals("CGT", S(e.getBreakendSequence()));
		assertEquals("", S(e.getAnchorSequence()));
		assertEquals(new BreakendSummary(0, BWD, 101, 100, 103), e.getBreakendSummary());
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
		realigned.setReadName(SplitReadIdentificationHelper.getSplitReadRealignments(r, false).get(0).getReadString());
		SplitReadIdentificationHelper.convertToSplitRead(r, ImmutableList.of(realigned));
		assertNull(r.getAttribute("SA"));
	}
	@Test
	public void getEvidenceID_paired_should_encode_read_direction_and_pair_info() {
		SAMRecord r = RP(1, 1, 100)[0];
		r.setReadName("ReadName");
		r.setCigarString("1S2M3S");
		r.setFirstOfPairFlag(true);
		SAMRecord r2 = r.deepCopy();
		r2.setFirstOfPairFlag(false);
		r2.setSecondOfPairFlag(true);
		assertEquals(4, Stream.of(
				SoftClipEvidence.create(SES(), BreakendDirection.Forward, r).getEvidenceID(),
				SoftClipEvidence.create(SES(), BreakendDirection.Backward, r).getEvidenceID(),
				SoftClipEvidence.create(SES(), BreakendDirection.Forward, r2).getEvidenceID(),
				SoftClipEvidence.create(SES(), BreakendDirection.Backward, r2).getEvidenceID()
				).distinct().count());
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
	public void getUntemplatedSequenceLength_should_match_soft_clip_with_no_realign() {
		SoftClipEvidence sc = SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "1S4M"));
		assertEquals(1, sc.getBreakendSequence().length);
		sc = SoftClipEvidence.create(SES(), BreakendDirection.Backward, Read(0, 1, "1S4M2S"));
		assertEquals(1, sc.getBreakendSequence().length);
		sc = SoftClipEvidence.create(SES(), BreakendDirection.Forward, Read(0, 1, "4M1S"));
		assertEquals(1, sc.getBreakendSequence().length);
		sc = SoftClipEvidence.create(SES(), BreakendDirection.Forward, Read(0, 1, "1S4M5S"));
		assertEquals(5, sc.getBreakendSequence().length);
	}
	@Test
	public void getLocalMapq_should_be_mapq() {
		assertEquals(5, SCE(FWD, withMapq(5, Read(0, 1, "4M2S"))).getLocalMapq());
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
		assertEquals(keep, !ses.shouldFilter(SoftClipEvidence.create(ses, FWD, r)));
		// reverse comp
		r = Read(0, 1000, String.format("%dS%dM", scLen, seq.length() - scLen));
		r.setReadBases(B(SequenceUtil.reverseComplement(seq)));
		r.setReadNegativeStrandFlag(true);
		assertEquals(keep, !ses.shouldFilter(SoftClipEvidence.create(ses, BWD, r)));
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
	public void isExact() {
		assertTrue(SCE(FWD, Read(0, 1, "1M1S")).isBreakendExact());
	}
	@Test
	public void getAnchorSequence_should_not_return_clipped_bases_on_other_side_of_the_read() {
		assertEquals("CGT", S(SCE(FWD, withSequence("ACGTNN", Read(0, 1, "1S3M2S"))[0]).getAnchorSequence()));
		assertEquals("CGT", S(SCE(BWD, withSequence("ACGTNN", Read(0, 1, "1S3M2S"))[0]).getAnchorSequence()));
	}
	@Test
	public void should_not_have_homology() {
		SoftClipEvidence e = SoftClipEvidence.create(SES(), FWD, withSequence("AAAAAAAAAA", Read(0, 10, "5M5S"))[0]);
		assertEquals("", e.getHomologySequence());
		assertEquals(0, e.getHomologyAnchoredBaseCount());
	}
}
