package au.edu.wehi.socrates;

import org.junit.Test;
import static org.junit.Assert.*;

import net.sf.samtools.SAMRecord;

public class NonReferenceReadPairTest extends TestHelper {
	public NonReferenceReadPair newPair(SAMRecord[] pair, int maxfragmentSize) {
		return new NonReferenceReadPair(pair[0], pair[1], maxfragmentSize);
	}
	@Test
	public void getLocalledMappedRead_should_return_local() {
		SAMRecord[] pair = OEA(0, 1, "100M", true);
		assertEquals(pair[0], newPair(pair, 1).getLocalledMappedRead());
	}
	@Test
	public void getRemoteReferenceIndex_should_return_remote() {
		SAMRecord[] pair = OEA(0, 1, "100M", true);
		assertEquals(pair[1], newPair(pair, 1).getNonReferenceRead());
	}
	@Test
	public void getBreakpointLocation_should_return_location_for_OEA() {
		assertTrue(newPair(OEA(0, 1, "100M", true), 1).getBreakpointLocation().getClass() == BreakpointLocation.class);
	}
	@Test
	public void getBreakpointLocation_should_return_interval_for_DP() {
		assertTrue(newPair(DP(0, 1, "100M", true, 0, 1, "100M", true), 1).getBreakpointLocation().getClass() == BreakpointInterval.class);
	}
	@Test
	public void getBreakpointLocation_forward_OEA_interval_should_allow_breakpoint_anywhere_in_fragment() {
		BreakpointLocation loc = newPair(OEA(0, 1, "2M3S", true), 10).getBreakpointLocation();
		// 123456789012345678901234567890
		// MMSSS
		// *--------* fragment
		//  |      |  possible breakpoint positions
		assertEquals(2, loc.start);
		assertEquals(9, loc.end);
		assertEquals(0, loc.referenceIndex);
	}
	@Test
	public void getBreakpointLocation_DP_calls_should_be_symmetrical() {
		SAMRecord[] pair = DP(0, 100, "1S3M1S", true, 1, 200, "5M", false); 
		BreakpointInterval loc1 = (BreakpointInterval)new NonReferenceReadPair(pair[0], pair[1], 20).getBreakpointLocation();
		BreakpointInterval loc2 = (BreakpointInterval)new NonReferenceReadPair(pair[1], pair[0], 20).getBreakpointLocation();
		assertEquals(loc1.qual, loc2.qual, 0);
		assertEquals(loc1.referenceIndex, loc2.referenceIndex2);
		assertEquals(loc1.start, loc2.start2);
		assertEquals(loc1.end, loc2.end2);
		assertEquals(loc1.direction, loc2.direction2);
		assertEquals(loc2.referenceIndex, loc1.referenceIndex2);
		assertEquals(loc2.start, loc1.start2);
		assertEquals(loc2.end, loc1.end2);
		assertEquals(loc2.direction, loc1.direction2);
	}
	@Test
	public void getBreakpointLocation_foward_OEA_interval_should_allow_breakpoint_anywhere_in_fragment() {
		BreakpointLocation loc = newPair(OEA(0, 1, "2M3S", true), 10).getBreakpointLocation();
		// 123456789012345678901234567890
		// MMSSS
		// *--------* fragment
		//  |      |  possible breakpoint positions
		assertEquals(2, loc.start);
		assertEquals(9, loc.end);
		assertEquals(0, loc.referenceIndex);
		assertEquals(BreakpointDirection.Forward, loc.direction);
	}
	@Test
	public void getBreakpointLocation_backward_OEA_interval_should_allow_breakpoint_anywhere_in_fragment() {
		BreakpointLocation loc = newPair(OEA(0, 10, "3S2M", false), 10).getBreakpointLocation();
		// 123456789012345678901234567890
		//       SSSMM
		//  *--------* fragment
		//   |      |   possible breakpoint positions
		assertEquals(3, loc.start);
		assertEquals(10, loc.end);
		assertEquals(0, loc.referenceIndex);
		assertEquals(BreakpointDirection.Backward, loc.direction);
	}
	private void dp_test_both(SAMRecord r1, SAMRecord r2, int maxFragmentSize,
			int expectedStart, int expectedEnd, int expectedReferenceIndex, BreakpointDirection expectedDirection) {
		BreakpointLocation loc = new NonReferenceReadPair(r1, r2, maxFragmentSize).getBreakpointLocation();
		assertEquals(expectedStart, loc.start);
		assertEquals(expectedEnd, loc.end);
		assertEquals(expectedReferenceIndex, loc.referenceIndex);
		assertEquals(expectedDirection, loc.direction);
		BreakpointInterval loc2 = (BreakpointInterval)new NonReferenceReadPair(r2, r1, maxFragmentSize).getBreakpointLocation();
		assertEquals(expectedStart, loc2.start2);
		assertEquals(expectedEnd, loc2.end2);
		assertEquals(expectedReferenceIndex, loc2.referenceIndex2);
		assertEquals(expectedDirection, loc2.direction2);
	}
	private void dp_test_both(SAMRecord[] dp, int maxFragmentSize,
			int expectedStart, int expectedEnd, int expectedReferenceIndex, BreakpointDirection expectedDirection) {
		dp_test_both(dp[0], dp[1], maxFragmentSize, expectedStart, expectedEnd, expectedReferenceIndex, expectedDirection);
	}
	@Test
	public void getBreakpointLocation_foward_backward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		// MMMMSS
		// *----------------------------*	max fragment
		//                         SMMMSS
		//    |                    | 		possible breakpoint positions
		dp_test_both(DP(0, 1, "4M2S", true, 1, 100, "1S3M2S", false), 30,
				4, 25, 0, BreakpointDirection.Forward);
	}
	@Test
	public void getBreakpointLocation_foward_forward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		// MMMMSS
		// *----------------------------*	max fragment
		//                         SSMMMS	need to flip the remote cigar since the breakpoint goes the other way
		//    |                     | 		possible breakpoint positions
		// [1, 31]							30
		// [6, 31] - read length			24
		// [4, 31] + local SC				26
		// [4, 25] - remote read length		20
		// [4, 27] + remote SC				22
		dp_test_both(DP(0, 1, "4M2S", true, 1, 100, "1S3M2S", true), 30,
				4, 26, 0, BreakpointDirection.Forward);
	}
	@Test
	public void getBreakpointLocation_backward_foward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		//                         SMMSSS
		// *----------------------------*	max fragment
		// SMMMSS
		//     |                    | 		possible breakpoint positions
		dp_test_both(DP(0, 26, "1S2M3S", false, 1, 100, "1S3M2S", true), 30,
				5, 26, 0, BreakpointDirection.Backward);
	}
	@Test
	public void getBreakpointLocation_backward_backward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		//                         SMMSSS
		// *----------------------------*	max fragment
		// SSMMMS							need to flip the remote cigar since the breakpoint goes the other way
		//      |                   | 		possible breakpoint positions
		dp_test_both(DP(0, 26, "1S2M3S", false, 1, 100, "1S3M2S", false), 30,
				6, 26, 0, BreakpointDirection.Backward);
	}
	@Test
	public void getBreakpointLocation_OEA_qual_should_match_mapq() {
		BreakpointLocation loc = newPair(OEA(0, 1, "2M3S", true), 10).getBreakpointLocation();
		assertEquals(10, loc.qual, 0);
	}

	@Test
	public void getBreakpointLocation_DP_qual_should_match_min_mapq() {
		SAMRecord[] pair = DP(0, 1, "100M", true, 0, 1, "100M", true);
		pair[0].setMappingQuality(1);
		pair[1].setMappingQuality(10);
		BreakpointInterval loc = (BreakpointInterval)newPair(pair, 1).getBreakpointLocation();
		assertEquals(1, loc.qual, 0);
		pair[0].setMappingQuality(10);
		pair[1].setMappingQuality(1);
		loc = (BreakpointInterval)newPair(pair, 1).getBreakpointLocation();
		assertEquals(1, loc.qual, 0);
	}
	public void getEvidenceID_should_match_read_name() {
		SAMRecord[] pair = DP(0, 1, "100M", true, 0, 1, "100M", true);
		pair[0].setReadName("EvidenceID");
		assertEquals("EvidenceID", newPair(pair, 1).getEvidenceID());
	}
}
