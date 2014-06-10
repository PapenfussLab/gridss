package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.EvidenceMetrics;


public class BreakpointSummaryTest extends TestHelper {
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
	@Test
	public void overlaps_should_consider_breakpoint_remote_location() {
		BreakpointSummary bp = new BreakpointSummary(0, FWD, 2, 3, 1, FWD, 2, 3, null); 
		for (BreakendSummary local : new BreakendSummary[] { 
				new BreakendSummary(0, FWD, 1, 4, null),
				new BreakendSummary(0, FWD, 2, 3, null),
				new BreakendSummary(0, FWD, 3, 4, null),
				new BreakendSummary(0, FWD, 1, 2, null)}) {
			for (BreakendSummary remote : new BreakendSummary[] { 
					new BreakendSummary(1, FWD, 1, 4, null),
					new BreakendSummary(1, FWD, 2, 3, null),
					new BreakendSummary(1, FWD, 3, 4, null),
					new BreakendSummary(1, FWD, 1, 2, null)}) {
				assertTrue(new BreakpointSummary(local, remote, null).overlaps(bp));
				assertTrue(bp.overlaps(new BreakpointSummary(local, remote, null)));
				// should not flip local and remote coordinates when checking for a match
				assertFalse(new BreakpointSummary(remote, local, null).overlaps(bp));
				assertFalse(bp.overlaps(new BreakpointSummary(remote, local, null)));
			}
		}
		// local does not overlap
		for (BreakendSummary local : new BreakendSummary[] { 
				new BreakendSummary(0, FWD, 1, 1, null),
				new BreakendSummary(0, FWD, 4, 4, null)}) {
			for (BreakendSummary remote : new BreakendSummary[] { 
					new BreakendSummary(1, FWD, 1, 4, null),
					new BreakendSummary(1, FWD, 2, 3, null),
					new BreakendSummary(1, FWD, 3, 4, null),
					new BreakendSummary(1, FWD, 1, 2, null),
					new BreakendSummary(1, FWD, 1, 1, null),
					new BreakendSummary(1, FWD, 4, 4, null)}) {
				assertFalse(new BreakpointSummary(local, remote, null).overlaps(bp));
				assertFalse(bp.overlaps(new BreakpointSummary(local, remote, null)));
			}
		}
		// remote does not overlap
		for (BreakendSummary local : new BreakendSummary[] { 
				new BreakendSummary(0, FWD, 1, 4, null),
				new BreakendSummary(0, FWD, 2, 3, null),
				new BreakendSummary(0, FWD, 3, 4, null),
				new BreakendSummary(0, FWD, 1, 2, null),
				new BreakendSummary(0, FWD, 1, 1, null),
				new BreakendSummary(0, FWD, 4, 4, null)}) {
			for (BreakendSummary remote : new BreakendSummary[] { 
					new BreakendSummary(1, FWD, 1, 1, null),
					new BreakendSummary(1, FWD, 4, 4, null)}) {
				assertFalse(new BreakpointSummary(local, remote, null).overlaps(bp));
				assertFalse(bp.overlaps(new BreakpointSummary(local, remote, null)));
			}
		}
	}
}
