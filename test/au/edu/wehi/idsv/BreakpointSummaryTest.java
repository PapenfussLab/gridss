package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.EvidenceMetrics;


public class BreakpointSummaryTest {
	@Test
	public void BreakpointSummary_explicit_constructor() {
		EvidenceMetrics em = new EvidenceMetrics();
		BreakpointSummary loc = new BreakpointSummary(1, BreakendDirection.Forward, 2, 3, 4, BreakendDirection.Backward, 5, 6, em);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
		assertSame(em, loc.evidence);
	}
	@Test
	public void BreakpointSummary_BreakendSummary_constructor() {
		EvidenceMetrics em1 = new EvidenceMetrics();
		EvidenceMetrics em2 = new EvidenceMetrics();
		EvidenceMetrics em3 = new EvidenceMetrics();
		BreakpointSummary loc = new BreakpointSummary(
				new BreakendSummary(1, BreakendDirection.Forward, 2, 3, em1),
				new BreakendSummary(4, BreakendDirection.Backward, 5, 6, em2),
				em3);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(4, loc.referenceIndex2);
		assertEquals(BreakendDirection.Backward, loc.direction2);
		assertEquals(5, loc.start2);
		assertEquals(6, loc.end2);
		assertSame(em3, loc.evidence);
	}
}
