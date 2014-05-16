package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BreakendSummaryTest {
	@Test
	public void BreakendSummary_Constructor() {
		EvidenceMetrics em = new EvidenceMetrics();
		BreakendSummary loc = new BreakendSummary(1, BreakendDirection.Forward, 2, 3, em);
		assertEquals(1, loc.referenceIndex);
		assertEquals(BreakendDirection.Forward, loc.direction);
		assertEquals(2, loc.start);
		assertEquals(3, loc.end);
		assertEquals(em, loc.evidence);
	}
}
