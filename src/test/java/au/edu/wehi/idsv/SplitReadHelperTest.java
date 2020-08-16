package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

public class SplitReadHelperTest extends TestHelper {
	@Test
	public void getSplitReadRealignments_should_split_start() {
		SAMRecord r = Read(0, 1, "2S4M");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		List<FastqRecord> result = SplitReadHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator(), (byte)0);
		
		assertEquals(1, result.size());
		assertEquals("AC", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("12")), result.get(0).getBaseQualityString());
	}
	@Test
	public void getSplitReadRealignments_should_split_end() {
		SAMRecord r = Read(0, 1, "3M3S");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		List<FastqRecord> result = SplitReadHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator(), (byte)0);
		
		assertEquals(1, result.size());
		assertEquals("TTC", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("456")), result.get(0).getBaseQualityString());
	}
	@Test
	public void getSplitReadRealignments_should_split_both_ends() {
		SAMRecord r = Read(0, 1, "2S1M3S");
		r.setReadBases(B("ACGTTC"));
		r.setBaseQualities(B("123456"));
		List<FastqRecord> result = SplitReadHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator(), (byte)0);
		
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
		List<FastqRecord> result = SplitReadHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator(), (byte)0);
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
		List<FastqRecord> result = SplitReadHelper.getSplitReadRealignments(r, false, getContext().getEvidenceIDGenerator(), (byte)0);
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
		List<FastqRecord> result = SplitReadHelper.getSplitReadRealignments(r, true, getContext().getEvidenceIDGenerator(), (byte)0);
		assertEquals(1, result.size());
		assertEquals("unique#4", result.get(0).getReadName());
		assertEquals("GT", result.get(0).getReadString());
		assertEquals(SAMUtils.phredToFastq(B("21")), result.get(0).getBaseQualityString());
	}
	@Test
	public void should_encode_alignment_unique_name_offset() {
		String encodedName = SplitReadHelper.getSplitAlignmentFastqName("r", 3);
		SAMRecord r2 = Read(0, 1, "3H2M2S");
		r2.setReadName(encodedName);
		assertEquals(3, SplitReadHelper.getRealignmentFirstAlignedBaseReadOffset(r2));
		assertEquals("r", SplitReadHelper.getOriginatingAlignmentUniqueName(r2));
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1), null, true);
		
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1), null, true);
		
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1), null, true);
		
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1), null, true);
		
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1), null, true);
		
		assertEquals("value", supp1.getStringAttribute("xx"));
		assertTrue(supp1.getReadPairedFlag());
		assertTrue(supp1.getProperPairFlag());
		assertFalse(supp1.getReadUnmappedFlag());
		assertTrue(supp1.getMateUnmappedFlag());
		assertFalse(supp1.getReadNegativeStrandFlag());
		assertTrue(supp1.getMateNegativeStrandFlag());
		assertTrue(supp1.getFirstOfPairFlag());
		assertTrue(supp1.getSecondOfPairFlag());
		assertTrue(supp1.isSecondaryAlignment());
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1), null, true);
		
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1, supp2), null, true);
		
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
		SplitReadHelper.convertToSplitRead(primary, ImmutableList.of(supp1, supp2), null, true);
		
		assertEquals(primary.getReadName(), supp1.getReadName());
		assertEquals(primary.getReadName(), supp2.getReadName());
	}
	@Test
	public void replaceAlignment_should_favour_overlapping_alignment() {
		SAMRecord r = Read(0, 1, "100M");
		SAMRecord r1 = withReadName(r.getReadName() + "#50", Read(1, 1, "50M"))[0];
		SAMRecord r2 = withReadName(r.getReadName() + "#0", Read(0, 1, "10M"))[0];
		SAMRecord result = SplitReadHelper.replaceAlignment(r, ImmutableList.of(r1, r2), true, true);
		Assert.assertEquals(result, r2);
		Assert.assertEquals("10M90S", r.getCigarString());
	}
	@Test
	public void replaceAlignment_should_process_supplmentary_alignments() {
		SAMRecord r = Read(0, 1, "100M");
		SAMRecord r1 = Read(1, 1, "50S50M");
		SAMRecord r2 = Read(0, 1, "10M90S");
		SAMRecord result = SplitReadHelper.replaceAlignment(r, ImmutableList.of(r1, r2), true, false);
		Assert.assertEquals(result, r2);
		Assert.assertEquals("10M90S", r.getCigarString());
	}
	@Test
	public void replaceAlignment_should_unmap_if_no_realignment_available() {
		SAMRecord r = Read(0, 1, "100M");
		SplitReadHelper.replaceAlignment(r, ImmutableList.of(), true, true);
		Assert.assertTrue(r.getReadUnmappedFlag());
	}
	@Test
	public void getEditDistanceDelta_should_compare_cigar_to_read() {
		Assert.assertArrayEquals(new int[] { 0, 1, 0, 1, 0, 0}, SplitReadHelper.getEditDistanceDelta(withSequence("ATATAA", Read(0, 100, "2S4M"))[0], getContext().getReference(), true));
		Assert.assertArrayEquals(new int[] { 0, 1, 1, 1, 0}, SplitReadHelper.getEditDistanceDelta(Read(0, 100, "1M3I1M"), getContext().getReference(), true));
	}
	@Test
	public void getEditDistanceDelta_should_include_hard_clipping() {
		Assert.assertArrayEquals(new int[] { 1,  1, 1, 0, 0, 0, 0, 1, 1}, SplitReadHelper.getEditDistanceDelta(Read(0, 100, "3H4M2H"), getContext().getReference(), true));
	}
	@Test
	public void getEditDistanceDelta_should_allocate_deletion_to_outermost_base() {
		Assert.assertArrayEquals(new int[] { 0, 2, 0}, SplitReadHelper.getEditDistanceDelta(Read(0, 100, "1M2D2M"), getContext().getReference(), true));
		Assert.assertArrayEquals(new int[] { 3, 1, 1}, SplitReadHelper.getEditDistanceDelta(withSequence("TTT", Read(0, 100, "1M2D2M"))[0], getContext().getReference(), false));
	}
	@Test
	public void getEditDistanceDelta_should_use_read_base_offsets() {
		Assert.assertArrayEquals(new int[] { 0, 2, 0}, SplitReadHelper.getEditDistanceDelta(onNegative(Read(0, 100, "1M2D2M"))[0], getContext().getReference(), false));
		Assert.assertArrayEquals(new int[] { 1, 1, 3}, SplitReadHelper.getEditDistanceDelta(onNegative(withSequence("TTT", Read(0, 100, "1M2D2M")))[0], getContext().getReference(), true));
	}
	@Test
	public void getEditDistanceDelta_should_consider_off_contig_as_mismatch() {
		Assert.assertArrayEquals(new int[] { 1, 1, 0}, SplitReadHelper.getEditDistanceDelta(Read(0, 1, "2S1M"), getContext().getReference(), false));
	}
	@Test
	public void adjustSplitLocationToMinimiseEditDistance_should_handle_overlapping_alignments() {
		// S====X=SSS
		// SSS=====SS
		SAMRecord left = withSequence("AAAAATAAAT", Read(0, 100, "1S6M3S"))[0];
		SAMRecord right = withSequence("TTAAAAAAAA", Read(0, 200, "3S5M2S"))[0];
		SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(left, right, getContext().getReference());
		Assert.assertEquals("1S4M5S", left.getCigarString());
		Assert.assertEquals(100, left.getAlignmentStart());
		Assert.assertEquals("5S3M2S", right.getCigarString());
		Assert.assertEquals(202, right.getAlignmentStart());
	}
	@Test
	public void adjustSplitLocationToMinimiseEditDistance_should_shift_right_to_best_position() {
		// S====X=SSS
		// SSS=====SS
		SAMRecord left = withSequence("AAATTT", Read(0, 100, "2M4S"))[0];
		SAMRecord right = withSequence("TTTAAA", Read(0, 200, "2S4M"))[0];
		SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(left, right, getContext().getReference());
		Assert.assertEquals("3M3S", left.getCigarString());
		Assert.assertEquals(100, left.getAlignmentStart());
		Assert.assertEquals("3S3M", right.getCigarString());
		Assert.assertEquals(201, right.getAlignmentStart());
	}
	@Test
	public void adjustSplitLocationToMinimiseEditDistance_should_shift_left_to_best_position() {
		// S====X=SSS
		// SSS=====SS
		SAMRecord left = withSequence("AAATTT", Read(0, 100, "5M1S"))[0];
		SAMRecord right = withSequence("TTTAAA", Read(0, 200, "1S5M"))[0];
		SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(left, right, getContext().getReference());
		Assert.assertEquals("3M3S", left.getCigarString());
		Assert.assertEquals("3S3M", right.getCigarString());
	}
	@Test
	public void should_adjustSplitLocationToMinimiseEditDistance_should_handle_negative_strand_alignment() {
		SAMRecord left = onNegative(withSequence("TTTAAA", Read(0, 100, "1S5M")))[0];
		SAMRecord right = withSequence("TTTAAA", Read(0, 200, "1S5M"))[0];
		SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(left, right, getContext().getReference());
		Assert.assertEquals("3S3M", left.getCigarString());
		Assert.assertEquals("3S3M", right.getCigarString());
	}
	@Test
	public void should_adjustSplitLocationToMinimiseEditDistance_should_not_unmap_read() {
		SAMRecord left = withSequence("CCCGGGCCC", Read(0, 100, "3S3M3S"))[0];
		SAMRecord right = withSequence("AAAAAAAAA", Read(0, 200, "3S6M"))[0];
		SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(left, right, getContext().getReference());
		Assert.assertEquals("3S1M5S", left.getCigarString());
		Assert.assertEquals("4S5M", right.getCigarString());
	}
	@Test
	public void adjustSplitLocationToMinimiseEditDistance_should_not_shift_optimal_alignments() {
		String seq = S(RANDOM).substring(0, 50) + S(RANDOM).substring(75, 100) + S(RANDOM).substring(150, 175);
		SAMRecord r1 = Read(2, 1, "50M50S");
		SAMRecord r2 = Read(2, 76, "50S25M25S");
		SAMRecord r3 = Read(2, 151, "75S25M");
		r1.setReadBases(B(seq));
		r2.setReadBases(B(seq));
		r3.setReadBases(B(seq));
		assertEquals(0, SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(r1, r2, SMALL_FA));
		assertEquals(0, SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(r2, r3, SMALL_FA));
		assertEquals(0, SplitReadHelper.adjustSplitLocationToMinimiseEditDistance(r1, r3, SMALL_FA));
	}
	@Test
	public void getEditDistanceDelta_should_force_mismatch_outside_of_chr_bounds() {
		SAMRecord r = withSequence("NNAAT", Read(0, 1, "2S3M"))[0];
		int[] distance = SplitReadHelper.getEditDistanceDelta(r, SMALL_FA, false);
		Assert.assertArrayEquals(new int[] {1, 1, 0, 0, 1}, distance);
	}
	@Test
	public void rewrite_anchor_should_updated_alignment() {
		SAMRecord r = withSequence("AAAAATTTTTTTTTTNNNNN", Read(0, 10, "5S10M5S"))[0];
		SAMRecord ra = onNegative(Read(1, 15, "10S4M6S"))[0];
		SplitReadHelper.rewriteAnchor(r, ra);
		assertEquals("10S4M6S", r.getCigarString());
		assertEquals(15, r.getAlignmentStart());
		assertEquals(1, (int)r.getReferenceIndex());
		assertEquals("NNNNNAAAAAAAAAATTTTT", S(r.getReadBases()));
		assertEquals("polyA,10,+,5S10M5S,10,0", r.getStringAttribute("OA"));
	}
	@Test
	public void should_truncate_anchor_alignment_to_not_extend_past_initial_anchor() {
		SAMRecord r = withSequence("AAAAATTTTTTTTTTNNNNN", Read(0, 10, "8S7M5S"))[0];
		SAMRecord ra = onNegative(Read(1, 15, "1S17M2S"))[0];
		SplitReadHelper.rewriteAnchor(r, ra);
		assertTrue(r.getReadNegativeStrandFlag());
		assertEquals("5S7M8S", r.getCigarString());
		assertEquals(19, r.getAlignmentStart());
	}
	@Test
	public void rewriteAnchor_should_prefer_overlapping() {
		SAMRecord r = Read(0, 100, "10S10M");
		SAMRecord ra1 = Read(0, 105, "10S1M9S");
		SAMRecord ra2 = Read(1, 15, "10S9M1S");
		ra1.setSupplementaryAlignmentFlag(true);
		SplitReadHelper.rewriteAnchor(r, ImmutableList.of(ra1, ra2));
		assertEquals("10S1M9S", r.getCigarString());
	}
	@Test
	public void rewriteAnchor_should_prefer_primary() {
		SAMRecord r = Read(0, 100, "10S10M");
		SAMRecord ra1 = Read(0, 105, "10S1M9S");
		SAMRecord ra2 = Read(0, 101, "10S9M1S");
		ra2.setSupplementaryAlignmentFlag(true);
		SplitReadHelper.rewriteAnchor(r, ImmutableList.of(ra1, ra2));
		assertEquals("10S1M9S", r.getCigarString());
	}
	@Test
	public void rewriteAnchor_should_prefer_longer_alignment() {
		SAMRecord r = Read(0, 100, "10S10M");
		SAMRecord ra1 = Read(0, 105, "10S1M9S");
		SAMRecord ra2 = Read(0, 101, "10S9M1S");
		ra2.setSupplementaryAlignmentFlag(true);
		ra1.setSupplementaryAlignmentFlag(true);
		SplitReadHelper.rewriteAnchor(r, ImmutableList.of(ra1, ra2));
		assertEquals("10S9M1S", r.getCigarString());
	}
}