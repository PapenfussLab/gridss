package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.Collections;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;


public class SequentialReferenceCoverageLookupTest extends TestHelper {
	public ReferenceCoverageLookup init(List<SAMRecord> reads, int windowSize) {
		Collections.sort(reads, new SAMRecordCoordinateComparator());
		return new SequentialReferenceCoverageLookup(reads.iterator(), IDSV(reads), new SAMFlagReadPairConcordanceCalculator(IDSV(reads)), windowSize);
	}
	@Test
	public void readsPairsSupportingNoBreakendAfter_should_return_first_read_in_pair() {
		List<SAMRecord> reads = L(
				RP(0, 10, 20, 5),
				RP(0, 11, 21, 5),
				RP(0, 12, 19, 5));
		ReferenceCoverageLookup lookup = init(reads, 1);
		for (int i = 1; i < 100; i++) {
			int r = lookup.readPairsSupportingNoBreakendAfter(0, i);
			//          1         2         3         4
			// 12345678901234567890123456789012345678901234567890
			//          AAAAA-----AAAAA
			//           AAAAA-----AAAAA
			//            AAAAA--AAAAA
			//              1233321
			assertEquals(
				i < 14 ? 0 :
				i <= 14 ? 1 :
				i <= 15 ? 2 :
				i <= 18 ? 3 :
				i <= 19 ? 2 :
				i <= 20 ? 1 :
				0, r);
		}
	}
	@Test
	public void overlapping_read_pairs_should_not_count_to_paired_coverage() {
		List<SAMRecord> reads = L(
				RP(0, 10, 10, 5),
				RP(0, 11, 15, 5),
				RP(0, 12, 17, 5));
		ReferenceCoverageLookup lookup = init(reads, 1);
		for (int i = 1; i < 100; i++) {
			int r = lookup.readPairsSupportingNoBreakendAfter(0, i);
			assertEquals(i == 16 ? 1 : 0, r);
		}
	}
	@Test
	public void should_consider_referenceIndex() {
		List<SAMRecord> reads = L(
				RP(0, 10, 20, 5),
				RP(1, 10, 20, 5),
				RP(1, 10, 20, 5));
		ReferenceCoverageLookup lookup = init(reads, 1);
		assertEquals(1, lookup.readsSupportingNoBreakendAfter(0, 10));
		assertEquals(1, lookup.readPairsSupportingNoBreakendAfter(0, 16));
		assertEquals(2, lookup.readsSupportingNoBreakendAfter(1, 10));
		assertEquals(2, lookup.readPairsSupportingNoBreakendAfter(1, 16));
	}
	@Test
	public void should_not_require_contiguous_calls() {
		List<SAMRecord> reads = L(
				RP(0, 10, 20, 5),
				RP(0, 11, 21, 5),
				RP(0, 12, 19, 5));
		ReferenceCoverageLookup lookup = init(reads, 1);
		assertEquals(3, lookup.readPairsSupportingNoBreakendAfter(0, 17)); 
		assertEquals(1, lookup.readPairsSupportingNoBreakendAfter(0, 20));
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_allow_out_of_order_position() {
		List<SAMRecord> reads = L(
				RP(0, 1, 75, 50),
				RP(0, 50, 150, 50));
		ReferenceCoverageLookup lookup = init(reads, 1);
		lookup.readPairsSupportingNoBreakendAfter(0, 55);
		lookup.readPairsSupportingNoBreakendAfter(0, 54);
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_not_allow_out_of_order_referenceIndex() {
		List<SAMRecord> reads = L(
				RP(0, 1, 50, 10),
				RP(1, 101, 150, 10));
		ReferenceCoverageLookup lookup = init(reads, 1);
		lookup.readPairsSupportingNoBreakendAfter(1, 100);
		lookup.readPairsSupportingNoBreakendAfter(0, 99);
	}
	@Test
	public void should_allow_out_of_order_of_window_size() {
		List<SAMRecord> reads = L(
				RP(0, 1, 50, 10),
				RP(0, 99, 150, 10));
		ReferenceCoverageLookup lookup = init(reads, 50);
		lookup.readsSupportingNoBreakendAfter(0, 100); 
		lookup.readsSupportingNoBreakendAfter(0, 51);
	}
	@Test
	public void readsSupportingNoBreakendAfter_should_return_reference_reads() {
		List<SAMRecord> reads = L(
				new SAMRecord[] {
				RP(0, 10, 50, 5)[0],
				Read(0, 10, "5M1S"),
				Read(0, 10, "1S5M")},
				OEA(0, 10, "5M", true),
				OEA(0, 10, "5M", false),
				DP(0, 10, "5M", true, 1, 10, "5M", false)
				);
		ReferenceCoverageLookup lookup = init(reads, 1);
		for (int i = 1; i < 50; i++) {
			int r = lookup.readsSupportingNoBreakendAfter(0, i);
			//          1         2         3         4
			// 12345678901234567890123456789012345678901234567890
			//          AAAAA     
			//          6666      
			assertEquals(i >= 10 && i < 14 ? 6 : 0, r);
		}
	}
	@Test
	public void should_consider_CIGAR() {
		List<SAMRecord> reads;
		ReferenceCoverageLookup lookup;
		// soft clips
		reads = L(Read(0, 2, "2S3M2S"));
		lookup = init(reads, 1);
		assertEquals(0, lookup.readPairsSupportingNoBreakendAfter(0, 1));
		//          1         2         3         4
		// 12345678901234567890123456789012345678901234567890
		//          AAAAA---------------
		//          SAAAS----------?????

		// Should handle indels in read
		//          1         2         3         4
		// 12345678901234567890123456789012345678901234567890
		//          AAAAA---------------
		//          AAIIIAAA----------?????
		//          ADA----------?????
		;
	}
	@Test
	public void readsPairsSupportingNoBreakendAfter_should_not_return_OEA() {
		List<SAMRecord> reads = L(
				OEA(0, 10, "5M", true),
				OEA(0, 10, "5M", false)
				);
		ReferenceCoverageLookup lookup = init(reads, 1);
		for (int i = 1; i < 100; i++) {
			assertEquals(0, lookup.readPairsSupportingNoBreakendAfter(0, i));
			assertEquals(0, lookup.readPairsSupportingNoBreakendAfter(0, i));
		}
	}
	@Test
	public void readsPairsSupportingNoBreakendAfter_should_not_return_Discordant() {
		ReferenceCoverageLookup lookup = init(L(
				DP(0, 10, "5M", true, 0, 18, "5M", true),
				DP(0, 18, "5M", true, 0, 10, "5M", true),
				DP(0, 10, "5M", false, 0, 18, "5M", true),
				DP(0, 10, "5M", true, 0, 18, "5M", false)), 1);
		for (int i = 1; i < 100; i++) {
			assertEquals(0, lookup.readPairsSupportingNoBreakendAfter(0, i));
		}
	}
}
