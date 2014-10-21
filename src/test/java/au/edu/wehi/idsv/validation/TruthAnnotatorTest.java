package au.edu.wehi.idsv.validation;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.TestHelper;


public class TruthAnnotatorTest extends TestHelper {
	@Test
	public void match_should_match_both_breakpoint_breakends() {
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			0, 0));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, BWD, 10, 10), new BreakendSummary(0, FWD, 1, 1)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				0, 0));
	}
	@Test
	public void match_should_match_local_breakend() {
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, FWD, 1, 1),
			0, 0));
	}
	@Test
	public void match_should_match_remote_breakend() {
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, BWD, 10, 10),
			0, 0));
	}
	@Test
	public void match_should_not_match_if_not_exact_match() {
		assertFalse(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 11, 11)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				0, 0));
	}
	@Test
	public void match_allow_any_interval_overlap() {
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 10), new BreakendSummary(0, BWD, 50, 60)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 5, 5), new BreakendSummary(0, BWD, 40, 70)),
				0, 0));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 10), new BreakendSummary(0, BWD, 50, 60)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 10, 11), new BreakendSummary(0, BWD, 40, 50)),
				0, 0));
	}
	/*
	 * If the bases we insert match the reference bases adjacent to the breakend, then we'll allow
	 * matches to those additional reference bases
	 */
	@Test
	public void match_should_allow_novel_base_insertions_to_match_reference() {
		// forward
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, FWD, 1, 1),
			4, 4));
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, FWD, 2, 2),
			4, 3));
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, FWD, 3, 3),
			4, 2));
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, FWD, 4, 4),
			4, 1));
		assertTrue(TruthAnnotator.matches(
			new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
			new BreakendSummary(0, FWD, 5, 5),
			4, 0));
		// backward
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakendSummary(0, BWD, 10, 10),
				4, 4));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakendSummary(0, BWD, 9, 9),
				4, 3));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakendSummary(0, BWD, 8, 8),
				4, 2));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakendSummary(0, BWD, 7, 7),
				4, 1));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakendSummary(0, BWD, 6, 6),
				4, 0));
		
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, BWD, 2, 2), new BreakendSummary(0, FWD, 10, 10)),
				new BreakendSummary(0, BWD, 1, 1),
				4, 3));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, BWD, 2, 2), new BreakendSummary(0, FWD, 10, 10)),
				new BreakendSummary(0, FWD, 11, 11),
				4, 3));
		
		// breakpoint should allow margin on either side
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				4, 4));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 2, 2), new BreakendSummary(0, BWD, 10, 10)),
				4, 3));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 9, 9)),
				4, 3));
		
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 6, 6)),
				4, 0));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 2, 2), new BreakendSummary(0, BWD, 7, 7)),
				4, 0));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 3, 3), new BreakendSummary(0, BWD, 8, 8)),
				4, 0));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 4, 4), new BreakendSummary(0, BWD, 9, 9)),
				4, 0));
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 10, 10)),
				new BreakpointSummary(new BreakendSummary(0, FWD, 5, 5), new BreakendSummary(0, BWD, 10, 10)),
				4, 0));
	}
	@Test
	public void match_should_limit_breakend_reference_homology_matching_to_16bp() {
		assertTrue(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 1000, 1000)),
				new BreakendSummary(0, FWD, 17, 17),
				17, 1));
		assertFalse(TruthAnnotator.matches(
				new BreakpointSummary(new BreakendSummary(0, FWD, 1, 1), new BreakendSummary(0, BWD, 1000, 1000)),
				new BreakendSummary(0, FWD, 18, 18),
				17, 0));
	}
	@Ignore
	@Test
	public void match_should_require_reference_and_untemplated_sequence_base_matching_for_breakend_reference_homology() {
		// TODO: implement this
		// low priority as only used for validation testing
	}
}
