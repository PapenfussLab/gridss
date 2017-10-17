package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

public class SplitReadIdentificationHelperTest extends TestHelper {
	@Test
	public void getSplitReadRealignments_should_split_start() {
		SAMRecord r = Read(0, 1, "2S4M");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		List<FastqRecord> result = SplitReadIdentificationHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator());
		
		assertEquals(1, result.size());
		assertEquals("AC", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("12")), result.get(0).getBaseQualityString());
	}
	@Test
	public void getSplitReadRealignments_should_split_end() {
		SAMRecord r = Read(0, 1, "3M3S");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		List<FastqRecord> result = SplitReadIdentificationHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator());
		
		assertEquals(1, result.size());
		assertEquals("TTC", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("456")), result.get(0).getBaseQualityString());
	}
	@Test
	public void getSplitReadRealignments_should_split_both_ends() {
		SAMRecord r = Read(0, 1, "2S1M3S");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		List<FastqRecord> result = SplitReadIdentificationHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator());
		
		assertEquals(2, result.size());
		assertEquals("AC", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("12")), result.get(0).getBaseQualityString());
		
		assertEquals("TTC", result.get(1).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("456")), result.get(1).getBaseQualityString());
	}
	@Test
	public void getSplitReadRealignments_name_should_encode_offset_from_start_of_read() {
		SAMRecord r = Read(0, 1, "2S1M3S");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		r.setReadName("r");
		List<FastqRecord> result = SplitReadIdentificationHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator());
		assertEquals(2, result.size());
		assertEquals(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(r) + "#0", result.get(0).getReadName());
		assertEquals(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(r) + "#3", result.get(1).getReadName());
	}
	@Test
	public void getSplitReadRealignments_should_consider_strand() {
		SAMRecord r = Read(0, 1, "2S1M3S");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		r.setReadName("r");
		r.setReadNegativeStrandFlag(true);
		List<FastqRecord> result = SplitReadIdentificationHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator());
		assertEquals(2, result.size());
		// SSSMSS cigar
		// GAACGT read
		// 654321 base qualities
		// 012345 offset
		assertEquals(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(r) + "#4", result.get(0).getReadName());
		assertEquals("GT", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("21")), result.get(0).getBaseQualityString());
		
		assertEquals(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(r) + "#0", result.get(1).getReadName());
		assertEquals("GAA", result.get(1).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("654")), result.get(1).getBaseQualityString());
	}
	@Test
	public void getSplitReadRealignments_should_chain_realignments() {
		// HHHMSS cigar
		// ???CGT read
		// ???321 base qualities
		// 012345 offset
		SAMRecord r = Read(0, 1, "1M2S");
		r.setReadName("unique#3");
		r.setReadBases(B("CGT"));
		r.setBaseQualities(B("321"));
		List<FastqRecord> result = SplitReadIdentificationHelper.getSplitReadRealignments(r, true, getContext().getEvidenceIDGenerator());
		assertEquals(1, result.size());
		assertEquals("unique#4", result.get(0).getReadName());
		assertEquals("GT", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("21")), result.get(0).getBaseQualityString());
	}
	@Test
	public void should_encode_alignment_unique_name_offset() {
		String encodedName = SplitReadIdentificationHelper.getSplitAlignmentFastqName("r", 3);
		SAMRecord r2 = Read(0, 1, "3H2M2S");
		r2.setReadName(encodedName);
		assertEquals(3, SplitReadIdentificationHelper.getRealignmentFirstAlignedBaseReadOffset(r2));
		assertEquals("r", SplitReadIdentificationHelper.getOriginatingAlignmentUniqueName(r2));
	}
	@Test
	public void convertToSplitRead_should_expand_clipping() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("CGTA"));
		primary.setBaseQualities(B("1234"));
		SAMRecord supp1 = Read(1, 2, "1M1S");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#1");
		supp1.setReadBases(B("GT"));
		supp1.setBaseQualities(B("23"));
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1));
		
		assertEquals("CGTA", S(supp1.getReadBases()));
		assertEquals("1234", S(supp1.getBaseQualities()));
		assertEquals("1S1M2S", supp1.getCigarString());
	}
	@Test
	public void convertToSplitRead_should_expand_clipping_with_negative_primary() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("CGTA"));
		primary.setBaseQualities(B("1234"));
		primary.setReadNegativeStrandFlag(true);
		// <---
		// CGTA
		// 1234
		// SSSM primary
		
		// --->
		// TACG
		// 4321
		// MSSS primary
		SAMRecord supp1 = Read(1, 2, "1M");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#2");
		supp1.setReadBases(B("C"));
		supp1.setBaseQualities(B("2"));
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1));
		
		assertEquals("TACG", S(supp1.getReadBases()));
		assertEquals("4321", S(supp1.getBaseQualities()));
		assertEquals("2S1M1S", supp1.getCigarString());
	}
	@Test
	public void convertToSplitRead_should_expand_clipping_with_negative_supp() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("CGTA"));
		primary.setBaseQualities(B("1234"));
		// --->
		// CGTA
		// 1234
		// SSSM primary
		// ??   supp
		
		// <---
		// TACG
		// 4321
		//   ??
		//   SM
		SAMRecord supp1 = Read(1, 2, "1S1M");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#0");
		supp1.setReadBases(B("CG"));
		supp1.setBaseQualities(B("21"));
		supp1.setReadNegativeStrandFlag(true);
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1));
		
		assertEquals("TACG", S(supp1.getReadBases()));
		assertEquals("4321", S(supp1.getBaseQualities()));
		assertEquals("3S1M", supp1.getCigarString());
	}
	@Test
	public void convertToSplitRead_should_copy_primary_attributes() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("CGTA"));
		primary.setBaseQualities(B("1234"));
		primary.setAttribute("xx", "value");
		SAMRecord supp1 = Read(1, 2, "1M1S");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#1");
		supp1.setReadBases(B("GT"));
		supp1.setBaseQualities(B("23"));
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1));
		
		assertEquals("value", supp1.getStringAttribute("xx"));
	}
	@Test
	public void convertToSplitRead_should_copy_primary_flags() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("CGTA"));
		primary.setBaseQualities(B("1234"));
		primary.setAttribute("xx", "value");
		primary.setFlags(-1);
		primary.setSupplementaryAlignmentFlag(false);
		primary.setReadUnmappedFlag(false);
		primary.setReadNegativeStrandFlag(true);
		SAMRecord supp1 = Read(1, 2, "1M1S");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#1");
		supp1.setReadBases(B("GT"));
		supp1.setBaseQualities(B("23"));
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1));
		
		assertEquals("value", supp1.getStringAttribute("xx"));
		assertTrue(supp1.getReadPairedFlag());
		assertTrue(supp1.getProperPairFlag());
		assertFalse(supp1.getReadUnmappedFlag());
		assertTrue(supp1.getMateUnmappedFlag());
		assertFalse(supp1.getReadNegativeStrandFlag());
		assertTrue(supp1.getMateNegativeStrandFlag());
		assertTrue(supp1.getFirstOfPairFlag());
		assertTrue(supp1.getSecondOfPairFlag());
		assertTrue(supp1.getNotPrimaryAlignmentFlag());
		assertTrue(supp1.getReadFailsVendorQualityCheckFlag());
		assertTrue(supp1.getDuplicateReadFlag());
		assertTrue(supp1.getSupplementaryAlignmentFlag());
	}
	@Test
	public void convertToSplitRead_should_copy_sam_fields() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("CGTA"));
		primary.setBaseQualities(B("1234"));
		primary.setMateReferenceIndex(2);
		primary.setMateAlignmentStart(7);
		SAMRecord supp1 = Read(1, 2, "1M1S");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#1");
		supp1.setReadBases(B("GT"));
		supp1.setBaseQualities(B("23"));
		supp1.setMappingQuality(17);
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1));
		
		assertEquals("r", supp1.getReadName());
		assertEquals(2, (int)supp1.getMateReferenceIndex());
		assertEquals(7, supp1.getMateAlignmentStart());
		assertEquals(17, supp1.getMappingQuality());
	}
	@Test
	public void convertToSplitRead_should_set_SA_tag_ordered_by_primary_then_read_offset() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("AAAA"));
		primary.setBaseQualities(B("1234"));
		primary.setMateReferenceIndex(2);
		primary.setMateAlignmentStart(7);
		primary.setMappingQuality(15);
		primary.setAttribute("NM", 8);
		
		SAMRecord supp1 = Read(0, 10, "1M1S");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#1");
		supp1.setReadBases(B("GT"));
		supp1.setBaseQualities(B("23"));
		supp1.setMappingQuality(17);
		supp1.setAttribute("NM", 9);
		
		SAMRecord supp2 = Read(0, 20, "1M");
		supp2.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#0");
		supp2.setReadBases(B("C"));
		supp2.setBaseQualities(B("1"));
		supp2.setMappingQuality(19);
		supp2.setAttribute("NM", 10);
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1, supp2));
		
		assertEquals("polyA,20,+,1M3S,19,10;polyA,10,+,1S1M2S,17,9", primary.getStringAttribute("SA"));
		assertEquals("polyA,1,+,3S1M,15,8;polyA,20,+,1M3S,19,10", supp1.getStringAttribute("SA"));
		assertEquals("polyA,1,+,3S1M,15,8;polyA,10,+,1S1M2S,17,9", supp2.getStringAttribute("SA"));
	}
	@Test
	public void supplementary_alignments_should_have_primary_read_name() {
		SAMRecord primary = Read(0, 1, "3S1M");
		primary.setReadName("r");
		primary.setReadBases(B("AAAA"));
		primary.setBaseQualities(B("1234"));
		primary.setMateReferenceIndex(2);
		primary.setMateAlignmentStart(7);
		primary.setMappingQuality(15);
		primary.setAttribute("NM", 8);
		
		SAMRecord supp1 = Read(0, 10, "1M1S");
		supp1.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#1");
		supp1.setReadBases(B("GT"));
		supp1.setBaseQualities(B("23"));
		supp1.setMappingQuality(17);
		supp1.setAttribute("NM", 9);
		
		SAMRecord supp2 = Read(0, 20, "1M");
		supp2.setReadName(getContext().getEvidenceIDGenerator().getAlignmentUniqueName(primary) + "#0");
		supp2.setReadBases(B("C"));
		supp2.setBaseQualities(B("1"));
		supp2.setMappingQuality(19);
		supp2.setAttribute("NM", 10);
		SplitReadIdentificationHelper.convertToSplitRead(primary, ImmutableList.of(supp1, supp2));
		
		assertEquals(primary.getReadName(), supp1.getReadName());
		assertEquals(primary.getReadName(), supp2.getReadName());
	}
}