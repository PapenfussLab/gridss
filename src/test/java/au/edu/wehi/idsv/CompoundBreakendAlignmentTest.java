package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import com.google.common.collect.ImmutableList;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class CompoundBreakendAlignmentTest extends TestHelper {
	@Test
	public void getSubsequentBreakpointAlignmentPairs_should_return_alignment_pairs_with_other_bases_soft_clipped() {
		// 1234567890
		// ATTCGTTTCATNAC
		// --SMMMMMSSSSSS
		//   MM
		//          SSMMM
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null,
				new BreakendSummary(0, FWD, 2, 2),
				B("AT"),
				B("12"),
				B("TCGTTTCATNAC"),
				B("234567890123"),
				ImmutableList.of(
						withName("0#0#0#readName", withMapq(40, withQual(B("234567890123"), withSequence(B("TCGTTTCATNAC"), Read(0, 10, "2S5M5S")))))[0],
						withName("0#0#0#readName", onNegative(withMapq(40, withQual(B("32"), withSequence(B("GA"), Read(1, 100, "2M"))))))[0],
						withName("0#0#7#readName", withMapq(40, withQual(B("90123"), withSequence(B("ATNAC"), Read(0, 200, "2S3M")))))[0]
				));
		List<Pair<SAMRecord, SAMRecord>> pairs = cba.getSubsequentBreakpointAlignmentPairs();
		assertEquals(2, pairs.size());
	}
	@Test
	public void should_apply_realignment_mapq_filter() {
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null, new BreakendSummary(0, FWD, 2, 2), B("AT"), B("12"), B("CGT"), B("123"), ImmutableList.of(
				withMapq(5, Read(0, 3, "3M"))[0]));
		assertTrue(cba.getSimpleBreakendRealignment().getReadUnmappedFlag());
		
		pc.getRealignmentParameters().mapqUniqueThreshold = 5;
		
		cba = new CompoundBreakendAlignment(pc, null, new BreakendSummary(0, FWD, 2, 2), B("AT"), B("12"), B("CGT"), B("123"), ImmutableList.of(
				withMapq(5, Read(0, 3, "3M"))[0]));
		assertFalse(cba.getSimpleBreakendRealignment().getReadUnmappedFlag());
	}
	@Test
	public void getSubsequentBreakpointAlignmentPairs_simple_FWD() {
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null,
				new BreakendSummary(0, FWD, 2, 2),
				B("AT"),
				B("12"),
				B("TCGTTTCATNACT"),
				B("2345678901234"),
				ImmutableList.of(
						withName("0#0#0#readName", withMapq(40, withQual(B("2345678901234"), withSequence(B("TCGTTTCATNACT"), Read(0, 10, "1S3M1D2M7S")))))[0],
						withName("0#0#8#readName", withMapq(40, withQual(B("01234"), withSequence(B("TNACT"), Read(0, 200, "1S4M")))))[0]
				));
		List<Pair<SAMRecord, SAMRecord>> pairs = cba.getSubsequentBreakpointAlignmentPairs();
		assertEquals(1, pairs.size());
		SAMRecord firstAnchor = pairs.get(0).getLeft();
		SAMRecord firstRealign = pairs.get(0).getRight();
		assertEquals(15, firstAnchor.getReadLength());
		assertEquals("ATTCGTTTCATNACT", S(firstAnchor.getReadBases()));
		assertEquals("122345678901234", S(firstAnchor.getBaseQualities()));
		assertEquals("3S3M1D2M7S", firstAnchor.getCigarString());
		assertEquals(10, firstAnchor.getAlignmentStart());
		
		assertEquals(7, firstRealign.getReadLength());
		assertEquals("CATNACT", S(firstRealign.getReadBases()));
		assertEquals("8901234", S(firstRealign.getBaseQualities()));
		assertEquals("3S4M", firstRealign.getCigarString());
		assertEquals(200, firstRealign.getAlignmentStart());
	}
	@Test
	public void getSubsequentBreakpointAlignmentPairs_single_bases_FWD() {
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null,
				new BreakendSummary(0, FWD, 2, 2),
				B("NN"),
				B("12"),
				B("ACGT"),
				B("1234"),
				ImmutableList.of(
						withName("0#0#0#readName", withMapq(40, withQual(B("A"), withSequence(B("1"), Read(0, 10, "1M")))))[0],
						withName("0#0#1#readName", withMapq(40, withQual(B("C"), withSequence(B("2"), Read(0, 10, "1M")))))[0],
						withName("0#0#2#readName", withMapq(40, withQual(B("G"), withSequence(B("3"), Read(0, 10, "1M")))))[0],
						withName("0#0#3#readName", withMapq(40, withQual(B("T"), withSequence(B("4"), Read(0, 10, "1M")))))[0]
				));
		List<Pair<SAMRecord, SAMRecord>> pairs = cba.getSubsequentBreakpointAlignmentPairs();
		assertEquals(3, pairs.size());
		assertEquals("NNACGT", S(pairs.get(0).getLeft().getReadBases()));
		assertEquals("2S1M3S", pairs.get(0).getLeft().getCigarString());
		assertEquals("CGT", S(pairs.get(0).getRight().getReadBases()));
		assertEquals("1M2S", pairs.get(0).getRight().getCigarString());
		
		assertEquals("NNACGT", S(pairs.get(1).getLeft().getReadBases()));
		assertEquals("3S1M2S", pairs.get(1).getLeft().getCigarString());
		assertEquals("GT", S(pairs.get(1).getRight().getReadBases()));
		assertEquals("1M1S", pairs.get(1).getRight().getCigarString());
		
		assertEquals("NNACGT", S(pairs.get(2).getLeft().getReadBases()));
		assertEquals("4S1M1S", pairs.get(2).getLeft().getCigarString());
		assertEquals("T", S(pairs.get(2).getRight().getReadBases()));
		assertEquals("1M", pairs.get(2).getRight().getCigarString());
	}
	@Test
	public void getSubsequentBreakpointAlignmentPairs_single_bases_BWD() {
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null,
				new BreakendSummary(0, BWD, 2, 2),
				B("NN"),
				B("12"),
				B("ACGT"),
				B("1234"),
				ImmutableList.of(
						withName("0#0#0#readName", withMapq(40, withQual(B("A"), withSequence(B("1"), Read(0, 10, "1M")))))[0],
						withName("0#0#1#readName", withMapq(40, withQual(B("C"), withSequence(B("2"), Read(0, 10, "1M")))))[0],
						withName("0#0#2#readName", withMapq(40, withQual(B("G"), withSequence(B("3"), Read(0, 10, "1M")))))[0],
						withName("0#0#3#readName", withMapq(40, withQual(B("T"), withSequence(B("4"), Read(0, 10, "1M")))))[0]
				));
		List<Pair<SAMRecord, SAMRecord>> pairs = cba.getSubsequentBreakpointAlignmentPairs();
		assertEquals(3, pairs.size());
		
		assertEquals("ACGTNN", S(pairs.get(0).getRight().getReadBases()));
		assertEquals("1S1M4S", pairs.get(0).getRight().getCigarString());
		assertEquals("A", S(pairs.get(0).getLeft().getReadBases()));
		assertEquals("1M", pairs.get(0).getLeft().getCigarString());
		
		assertEquals("ACGTNN", S(pairs.get(1).getRight().getReadBases()));
		assertEquals("2S1M3S", pairs.get(1).getRight().getCigarString());
		assertEquals("AC", S(pairs.get(1).getLeft().getReadBases()));
		assertEquals("1S1M", pairs.get(1).getLeft().getCigarString());
		
		assertEquals("ACGTNN", S(pairs.get(2).getRight().getReadBases()));
		assertEquals("3S1M2S", pairs.get(2).getRight().getCigarString());
		assertEquals("ACG", S(pairs.get(2).getLeft().getReadBases()));
		assertEquals("2S1M", pairs.get(2).getLeft().getCigarString());
	}
	@Test
	public void getSubsequentBreakpointAlignmentPairs_overlapping_single_bases_FWD() {
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null,
				new BreakendSummary(0, FWD, 2, 2),
				B("NN"),
				B("12"),
				B("ACGT"),
				B("1234"),
				ImmutableList.of(
						withName("0#0#0#readName", withMapq(40, withQual(B("ACGT"), withSequence(B("1234"), Read(0, 10, "1M3S")))))[0],
						withName("0#0#0#readName", withMapq(40, withQual(B("ACGT"), withSequence(B("1234"), Read(0, 10, "1S1M2S")))))[0],
						withName("0#0#0#readName", withMapq(40, withQual(B("ACGT"), withSequence(B("1234"), Read(0, 10, "2S1M1S")))))[0],
						withName("0#0#0#readName", withMapq(40, withQual(B("ACGT"), withSequence(B("1234"), Read(0, 10, "3S1M")))))[0]
				));
		List<Pair<SAMRecord, SAMRecord>> pairs = cba.getSubsequentBreakpointAlignmentPairs();
		assertEquals(3, pairs.size());
		assertEquals("NNACGT", S(pairs.get(0).getLeft().getReadBases()));
		assertEquals("2S1M3S", pairs.get(0).getLeft().getCigarString());
		assertEquals("CGT", S(pairs.get(0).getRight().getReadBases()));
		assertEquals("1M2S", pairs.get(0).getRight().getCigarString());
		
		assertEquals("NNACGT", S(pairs.get(1).getLeft().getReadBases()));
		assertEquals("3S1M2S", pairs.get(1).getLeft().getCigarString());
		assertEquals("GT", S(pairs.get(1).getRight().getReadBases()));
		assertEquals("1M1S", pairs.get(1).getRight().getCigarString());
		
		assertEquals("NNACGT", S(pairs.get(2).getLeft().getReadBases()));
		assertEquals("4S1M1S", pairs.get(2).getLeft().getCigarString());
		assertEquals("T", S(pairs.get(2).getRight().getReadBases()));
		assertEquals("1M", pairs.get(2).getRight().getCigarString());
	}
	@Test
	public void simple_compound_test_fwd() {
		//          1         2         3         4         5         6         7      
		// 123456789012345678901234567890123456789012345678901234567890123456789012345
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA
		// *****SSSSSSSSSSSSSSSSSSSSSSSSSMMMMMSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
		//          SSSSSSSSMMMMMMMMSSSSSSSSSSSSSSSSS
		//                                        M       
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(getContext(), null,
				new BreakendSummary(0, FWD, 5, 5),
				B("CATTA"),
				B("00000"),
				B("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA"),
				B("1111111111111111111111111111111111111111111111111111111111111111111111"),
				ImmutableList.of(
					withReadName("0#0#0#readname", withSequence(B("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA"),
							withQual(B(40,"ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA".length()), Read(1, 100, "25S5M40S"))))[0],
					withReadName("0#0#4#readname", withSequence(B("CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT"),
							withQual(B(40,"CAAGAGCGGGTTGTATTCGACGCCAAGTCAGCT".length()), Read(2, 200, "8S8M17S"))))[0],
					withReadName("0#0#34#readname", withQual(B("1"), withSequence("G", Read(0, 40, "1M"))))[0]
					));
		assertEquals("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", S(cba.getSimpleBreakendRealignment().getReadBases()));
		assertEquals("12S8M50S", cba.getSimpleBreakendRealignment().getCigarString());
		assertEquals(200, cba.getSimpleBreakendRealignment().getAlignmentStart());
		assertEquals(2, (int)cba.getSimpleBreakendRealignment().getReferenceIndex());
		assertEquals(2, cba.getSubsequentBreakpointAlignmentPairs().size());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", S(cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getReadBases()));
		assertEquals("TCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", S(cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getReadBases()));
		assertEquals(2, (int)cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getReferenceIndex());
		assertEquals(200, cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getAlignmentStart());
		assertEquals("17S8M50S", cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getCigarString());
		assertFalse(cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getReadNegativeStrandFlag());
		assertEquals(1, (int)cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getReferenceIndex());
		assertEquals(100, cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getAlignmentStart());
		assertEquals("5S5M40S", cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getCigarString());
		assertFalse(cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getReadNegativeStrandFlag());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", S(cba.getSubsequentBreakpointAlignmentPairs().get(1).getLeft().getReadBases()));
		assertEquals("GTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAA", S(cba.getSubsequentBreakpointAlignmentPairs().get(1).getRight().getReadBases()));
		assertEquals(1, (int)cba.getSubsequentBreakpointAlignmentPairs().get(1).getLeft().getReferenceIndex());
		assertEquals(100, cba.getSubsequentBreakpointAlignmentPairs().get(1).getLeft().getAlignmentStart());
		assertEquals("30S5M40S", cba.getSubsequentBreakpointAlignmentPairs().get(1).getLeft().getCigarString());
		assertFalse(cba.getSubsequentBreakpointAlignmentPairs().get(1).getRight().getReadNegativeStrandFlag());
		assertEquals(0, (int)cba.getSubsequentBreakpointAlignmentPairs().get(1).getRight().getReferenceIndex());
		assertEquals(40, cba.getSubsequentBreakpointAlignmentPairs().get(1).getRight().getAlignmentStart());
		assertEquals("4S1M35S", cba.getSubsequentBreakpointAlignmentPairs().get(1).getRight().getCigarString());
		assertFalse(cba.getSubsequentBreakpointAlignmentPairs().get(1).getRight().getReadNegativeStrandFlag());
	}
	@Test
	public void compound_test_fwd_fwd_bwd() {
		//          1         2         3         4      
		// 1234567890123456789012345678901234567890
		// CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAG
		// *****    SSSSMMMMMSSSSSS
		//                SSSSSSSMMMMMMMMSSSSSSSSS (but reversed)
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(getContext(), null,
				new BreakendSummary(0, FWD, 5, 5),
				B("CATTA"),
				B("00000"),
				B("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAG"),
				B("35353535353535353535353535353535353"),
				ImmutableList.of(
					withReadName("0#0#4#readname", withSequence(B("CAAGAGCGGGTTGTA"),
							withQual(B(40,"CAAGAGCGGGTTGTA".length()), Read(1, 100, "4S5M6S"))))[0],
					onNegative(withReadName("0#0#10#readname", withSequence(B(SequenceUtil.reverseComplement("CGGGTTGTATTCGACGCCAAGTCA")),
							withQual(B(40,"CGGGTTGTATTCGACGCCAAGTCA".length()), Read(2, 200, "9S8M7S")))))[0]
					));
		assertEquals("ATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAG", S(cba.getSimpleBreakendRealignment().getReadBases()));
		assertEquals("8S5M22S", cba.getSimpleBreakendRealignment().getCigarString());
		assertEquals(1, (int)cba.getSimpleBreakendRealignment().getReferenceIndex());
		assertEquals(100, cba.getSimpleBreakendRealignment().getAlignmentStart());
		assertEquals(1, cba.getSubsequentBreakpointAlignmentPairs().size());
		
		assertEquals(1, (int)cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getReferenceIndex());
		assertEquals(100, cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getAlignmentStart());
		assertEquals("13S5M22S", cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getCigarString());
		assertEquals("CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAG", S(cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getReadBases()));
		assertFalse(cba.getSubsequentBreakpointAlignmentPairs().get(0).getLeft().getReadNegativeStrandFlag());
		
		assertEquals(2, (int)cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getReferenceIndex());
		assertEquals(200, cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getAlignmentStart());
		assertEquals("10S8M4S", cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getCigarString());
		assertEquals(SequenceUtil.reverseComplement("GTTGTATTCGACGCCAAGTCAG"), S(cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getReadBases()));
		assertTrue(cba.getSubsequentBreakpointAlignmentPairs().get(0).getRight().getReadNegativeStrandFlag());
	}
}
