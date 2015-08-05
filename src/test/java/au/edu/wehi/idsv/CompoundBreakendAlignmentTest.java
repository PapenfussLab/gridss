package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import com.google.common.collect.ImmutableList;

import htsjdk.samtools.SAMRecord;

public class CompoundBreakendAlignmentTest extends TestHelper {
	@Test
	public void getSubsequentBreakpointAlignmentPairs_should_return_alignment_pairs_with_other_bases_soft_clipped() {
		ProcessingContext pc = getContext();
		CompoundBreakendAlignment cba = new CompoundBreakendAlignment(pc, null,
				new BreakendSummary(0, FWD, 2, 2),
				B("AT"),
				B("12"),
				B("TCGTTTCATNAC"),
				B("234567890123"),
				ImmutableList.of(
						withName("0#0#0#readName", withMapq(40, withQual(B("234567890123"), withSequence(B("TCGTTTCATNAC"), Read(0, 10, "1S5M5S")))))[0],
						withName("0#0#0#readName", onNegative(withMapq(40, withQual(B("32"), withSequence(B("GA"), Read(1, 100, "2M"))))))[0],
						withName("0#0#8#readName", withMapq(40, withQual(B("90123"), withSequence(B("ATNAC"), Read(0, 200, "2S3M")))))[0]
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
}
