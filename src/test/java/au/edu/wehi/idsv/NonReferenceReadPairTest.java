package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

public class NonReferenceReadPairTest extends TestHelper {
	public NonReferenceReadPair newPair(SAMRecord[] pair, int maxfragmentSize) {
		return NonReferenceReadPair.create(pair[0], pair[1], SES(maxfragmentSize));
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_abort_if_max_frag_size_not_sane() {
		SAMRecord[] pair = OEA(0, 1, "100M", true);
		newPair(pair, 99);
	}
	@Test
	public void getLocalledMappedRead_should_return_local() {
		SAMRecord[] pair = OEA(0, 1, "100M", true);
		assertEquals(pair[0], newPair(pair, 100).getLocalledMappedRead());
	}
	@Test
	public void getRemoteReferenceIndex_should_return_remote() {
		SAMRecord[] pair = OEA(0, 1, "100M", true);
		assertEquals(pair[1], newPair(pair, 100).getNonReferenceRead());
	}
	@Test
	public void getBreakendSummary_should_return_location_for_OEA() {
		assertTrue(newPair(OEA(0, 1, "100M", true), 300).getBreakendSummary().getClass() == BreakendSummary.class);
	}
	@Test
	public void getBreakendSummary_should_return_interval_for_DP() {
		assertTrue(newPair(DP(0, 1, "100M", true, 0, 100, "100M", true), 300).getBreakendSummary().getClass() == BreakpointSummary.class);
	}
	@Test
	public void getBreakendSummary_forward_OEA_interval_should_allow_breakpoint_anywhere_in_fragment() {
		BreakendSummary loc = newPair(OEA(0, 1, "2M3S", true), 10).getBreakendSummary();
		// 123456789012345678901234567890
		// MMSSS
		// *--------* fragment
		//  |      |  possible breakpoint positions
		assertEquals(2, loc.start);
		assertEquals(9, loc.end);
		assertEquals(0, loc.referenceIndex);
	}
	@Test
	public void getBreakendSummary_DP_calls_should_be_symmetrical() {
		SAMRecord[] pair = DP(0, 100, "1S3M1S", true, 1, 200, "5M", false); 
		BreakpointSummary loc1 = (BreakpointSummary)NonReferenceReadPair.create(pair[0], pair[1], SES(20)).getBreakendSummary();
		BreakpointSummary loc2 = (BreakpointSummary)NonReferenceReadPair.create(pair[1], pair[0], SES(20)).getBreakendSummary();
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
	public void getBreakendSummary_foward_OEA_interval_should_allow_breakpoint_anywhere_in_fragment() {
		BreakendSummary loc = newPair(OEA(0, 1, "2M3S", true), 10).getBreakendSummary();
		// 123456789012345678901234567890
		// MMSSS
		// *--------* fragment
		//  |      |  possible breakpoint positions
		assertEquals(2, loc.start);
		assertEquals(9, loc.end);
		assertEquals(0, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
	}
	@Test
	public void getBreakendSummary_backward_OEA_interval_should_allow_breakpoint_anywhere_in_fragment() {
		BreakendSummary loc = newPair(OEA(0, 10, "3S2M", false), 10).getBreakendSummary();
		// 123456789012345678901234567890
		//       SSSMM
		//  *--------* fragment
		//   |      |   possible breakpoint positions
		assertEquals(3, loc.start);
		assertEquals(10, loc.end);
		assertEquals(0, loc.referenceIndex);
		assertEquals(BreakendDirection.Backward, loc.direction);
	}
	private void dp_test_both(SAMRecord r1, SAMRecord r2, int maxFragmentSize,
			int expectedStart, int expectedEnd, int expectedReferenceIndex, BreakendDirection expectedDirection) {
		BreakendSummary loc = NonReferenceReadPair.create(r1, r2, SES(maxFragmentSize)).getBreakendSummary();
		assertEquals(expectedStart, loc.start);
		assertEquals(expectedEnd, loc.end);
		assertEquals(expectedReferenceIndex, loc.referenceIndex);
		assertEquals(expectedDirection, loc.direction);
		BreakpointSummary loc2 = (BreakpointSummary)NonReferenceReadPair.create(r2, r1, SES(maxFragmentSize)).getBreakendSummary();
		assertEquals(expectedStart, loc2.start2);
		assertEquals(expectedEnd, loc2.end2);
		assertEquals(expectedReferenceIndex, loc2.referenceIndex2);
		assertEquals(expectedDirection, loc2.direction2);
	}
	private void dp_test_both(SAMRecord[] dp, int maxFragmentSize,
			int expectedStart, int expectedEnd, int expectedReferenceIndex, BreakendDirection expectedDirection) {
		dp_test_both(dp[0], dp[1], maxFragmentSize, expectedStart, expectedEnd, expectedReferenceIndex, expectedDirection);
	}
	@Test
	public void getBreakendSummary_foward_backward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		// MMMMSS
		// *----------------------------*	max fragment
		//                         SMMMSS
		//    |                    | 		possible breakpoint positions
		dp_test_both(DP(0, 1, "4M2S", true, 1, 100, "1S3M2S", false), 30,
				4, 25, 0, BreakendDirection.Forward);
	}
	@Test
	public void getBreakendSummary_foward_forward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
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
				4, 26, 0, BreakendDirection.Forward);
	}
	@Test
	public void getBreakendSummary_backward_foward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		//                         SMMSSS
		// *----------------------------*	max fragment
		// SMMMSS
		//     |                    | 		possible breakpoint positions
		dp_test_both(DP(0, 26, "1S2M3S", false, 1, 100, "1S3M2S", true), 30,
				5, 26, 0, BreakendDirection.Backward);
	}
	@Test
	public void getBreakendSummary_backward_backward_DP_interval_should_allow_breakpoint_in_unsequenced_content() {
		//          1         2         3         4         5         6
		// 123456789012345678901234567890123456789012345678901234567890
		//                         SMMSSS
		// *----------------------------*	max fragment
		// SSMMMS							need to flip the remote cigar since the breakpoint goes the other way
		//      |                   | 		possible breakpoint positions
		dp_test_both(DP(0, 26, "1S2M3S", false, 1, 100, "1S3M2S", false), 30,
				6, 26, 0, BreakendDirection.Backward);
	}
	@Test
	public void getBreakendSummary_OEA_qual_should_match_mapq() {
		BreakendSummary loc = newPair(OEA(0, 1, "2M3S", true), 10).getBreakendSummary();
	}
	@Test
	public void DP_should_set_dp_evidence() {
		SAMRecord[] pair = DP(0, 1, "100M", true, 0, 1, "100M", true);
		pair[0].setMappingQuality(1);
		pair[1].setMappingQuality(10);
		BreakpointSummary loc = (BreakpointSummary)newPair(pair, 300).getBreakendSummary();
		pair[0].setMappingQuality(10);
		pair[1].setMappingQuality(1);
		loc = (BreakpointSummary)newPair(pair, 300).getBreakendSummary();
	}
	public void getEvidenceID_should_match_read_name() {
		SAMRecord[] pair = DP(0, 1, "100M", true, 0, 1, "100M", true);
		pair[0].setReadName("EvidenceID");
		assertEquals("EvidenceID", newPair(pair, 1).getEvidenceID());
	}
	@Test
	public void getBreakendSummary_DP_should_restrict_to_within_fragment_if_concordant_placement_but_not_length() {
		NonReferenceReadPair pair = newPair(DP(0, 101, "3M", true, 0, 110, "3M", false), 100);
		// 1
		// 0
		// 012345678901234567890
		//  ==>      <== only support a BP within this interval
		//    ^------ forward
		//     ------^ backward
		assertEquals(103, pair.getBreakendSummary().start);
		assertEquals(109, pair.getBreakendSummary().end);
		assertEquals(104, ((BreakpointSummary)pair.getBreakendSummary()).start2);
		assertEquals(110, ((BreakpointSummary)pair.getBreakendSummary()).end2);
	}
	@Test
	public void getBreakendSummary_DP_should_not_restrict_on_alt_contigs() {
		NonReferenceReadPair pair = newPair(DP(0, 101, "3M", true, 1, 110, "3M", false), 100);
		assertEquals(103, pair.getBreakendSummary().start);
		assertEquals(103 + (100 - 3 - 3), pair.getBreakendSummary().end);
		assertEquals(110, ((BreakpointSummary)pair.getBreakendSummary()).end2);
		assertEquals(110 - (100 - 3 - 3), ((BreakpointSummary)pair.getBreakendSummary()).start2);
	}
	@Test
	public void overlapping_PE_reads_should_provide_no_breakpoint_support() {
		// 012345678901234567890
		//  ==>
		//    <==
		assertNull(newPair(DP(0, 1, "3M", true, 0, 3, "3M", false), 100).getBreakendSummary());
		assertNull(newPair(DP(0, 3, "3M", false, 0, 1, "3M", true), 100).getBreakendSummary());
	}
	@Test
	public void overlapping_PE_reads_should_suppress_tandem_duplication_support() {
		// if we map over the ends of each other, we're not just a short fragment
		// 012345678901234567890
		//     =====>
		//  <=====
		assertNotNull(newPair(DP(0, 4, "6M", true, 0, 1, "6M", false), 100).getBreakendSummary());
		assertNotNull(newPair(DP(0, 1, "6M", false, 0, 4, "6M", true), 100).getBreakendSummary());
	}
	@Test
	public void should_handle_concordant_position_but_discordant_chr() {
		newPair(DP(1, 1, "100M", true, 2, 5, "100M", true), 300);
		newPair(DP(1, 2, "100M", true, 2, 4, "100M", true), 300);
		newPair(DP(1, 3, "100M", true, 2, 6, "100M", true), 300);
		newPair(OEA(1, 4, "100M", false), 300);
	}
	@Test
	public void create_should_instanciate_derived_classes() {
		assertTrue(newPair(DP(1, 1, "100M", true, 2, 5, "100M", true), 300) instanceof DiscordantReadPair);
		assertTrue(newPair(OEA(1, 4, "100M", false), 300) instanceof UnmappedMateReadPair);
	}
	@Test
	public void getLocalMapq_should_be_anchored_mapq() {
		SAMRecord[] pair = DP(1, 1, "100M", true, 2, 5, "100M", true);
		pair[0].setMappingQuality(5);
		pair[1].setMappingQuality(10);
		assertEquals(5, newPair(pair, 300).getLocalMapq());
	}
	@Test
	public void getLocalBaseLength_should_be_read_length() {
		assertEquals(100, newPair(OEA(0, 1, "100M", true), 300).getLocalBaseLength());
	}
	@Test
	public void getLocalBaseCount_should_be_read_length() {
		assertEquals(100, newPair(OEA(0, 1, "100M", true), 300).getLocalBaseCount());
	}
	@Test
	public void getLocalMaxBaseQual_local_mapped_quals() {
		SAMRecord[] pair = DP(1, 1, "3M1S", true, 2, 5, "1S3M", true);
		withQual(new byte[] { 1, 2, 3, 4}, pair[0]);
		withQual(new byte[] { 4, 5, 6, 7}, pair[1]);
		assertEquals(3, newPair(pair, 300).getLocalMaxBaseQual());
	}
	@Test
	public void getLocalTotalBaseQual_local_mapped_quals() {
		SAMRecord[] pair = DP(1, 1, "3M1S", true, 2, 5, "1S3M", true);
		withQual(new byte[] { 1, 2, 3, 4}, pair[0]);
		withQual(new byte[] { 4, 5, 6, 7}, pair[1]);
		assertEquals(1+2+3, newPair(pair, 300).getLocalTotalBaseQual());
	}
	@Test
	public void getBreakendSequence_should_be_null() {
		assertNull(newPair(OEA(0, 1, "100M", true), 300).getBreakendSequence());
	}
	@Test
	public void getBreakendQuality_should_be_null() {
		assertNull(newPair(OEA(0, 1, "100M", true), 300).getBreakendSequence());
	}
}
