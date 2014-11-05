package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;


public class EvidenceClusterSubsetProcessorTest extends TestHelper {
	public class TestEvidenceClusterSubsetProcessor extends EvidenceClusterSubsetProcessor {
		public TestEvidenceClusterSubsetProcessor(int fromReferenceIndex, int toReferenceIndex) {
				super(getContext(), fromReferenceIndex, toReferenceIndex);
		}
		@Override
		public boolean filterOut(BreakendSummary loc) {
			return super.filterOut(loc);
		}
		@Override
		public boolean filterOut(BreakpointSummary loc) {
			return super.filterOut(loc);
		}
	}
	@Test
	public void should_filterOut_breakend_on_alt_chr() {
		assertFalse(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakendSummary(0, FWD, 1, 1)));
		assertFalse(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakendSummary(1, FWD, 1, 1)));
		assertTrue(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakendSummary(2, FWD, 1, 1)));
	}
	@Test
	public void should_filterOut_breakpoint_not_between_chr_pair() {
		assertFalse(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(0,  FWD,  1,  1,  1,  FWD,  1,  1)));
		assertFalse(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(1,  FWD,  1,  1,  0,  FWD,  1,  1)));
		assertTrue(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(0,  FWD,  1,  1,  0,  FWD,  1,  1)));
		assertTrue(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(1,  FWD,  1,  1,  1,  FWD,  1,  1)));
		assertTrue(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(0,  FWD,  1,  1,  2,  FWD,  1,  1)));
		assertTrue(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(2,  FWD,  1,  1,  0,  FWD,  1,  1)));
		assertTrue(new TestEvidenceClusterSubsetProcessor(0, 1).filterOut(new BreakpointSummary(2,  FWD,  1,  1,  1,  FWD,  1,  1)));
	}
}
