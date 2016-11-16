package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.Collections;

import org.junit.Test;


public class EvidenceClusterSubsetProcessorTest extends TestHelper {
	public class TestEvidenceClusterSubsetProcessor extends EvidenceClusterSubsetProcessor {
		private ChromosomeFilteringIterator testIt;
		public TestEvidenceClusterSubsetProcessor(int fromReferenceIndex, int toReferenceIndex) {
			super(getContext(), Collections.<DirectedEvidence>emptyIterator(), fromReferenceIndex, toReferenceIndex);
			this.testIt = new ChromosomeFilteringIterator(Collections.<DirectedEvidence>emptyIterator(), fromReferenceIndex, toReferenceIndex);
		}
	}
	public void go(BreakendSummary be, boolean expectFiltered) {
		TestEvidenceClusterSubsetProcessor p = new TestEvidenceClusterSubsetProcessor(0, 1);
		assertEquals(expectFiltered, p.testIt.isFiltered(be));
		p.close();
	}
	public void go(BreakpointSummary bp, boolean expectFiltered) {
		TestEvidenceClusterSubsetProcessor p = new TestEvidenceClusterSubsetProcessor(0, 1);
		assertEquals(expectFiltered, p.testIt.isFiltered(bp));
		p.close();
	}
	@Test
	public void should_filterOut_breakend_on_alt_chr() {
		go(new BreakendSummary(0, FWD, 1), false);
		go(new BreakendSummary(1, FWD, 1), false);
		go(new BreakendSummary(2, FWD, 1), true);
	}
	@Test
	public void should_filterOut_breakpoint_not_between_chr_pair() {
		go(new BreakpointSummary(0,  FWD,  1,  1,  FWD,  1), false);
		go(new BreakpointSummary(1,  FWD,  1,  0,  FWD,  1), false);
		go(new BreakpointSummary(0,  FWD,  1,  0,  FWD,  1), true);
		go(new BreakpointSummary(1,  FWD,  1,  1,  FWD,  1), true);
		go(new BreakpointSummary(0,  FWD,  1,  2,  FWD,  1), true);
		go(new BreakpointSummary(2,  FWD,  1,  0,  FWD,  1), true);
		go(new BreakpointSummary(2,  FWD,  1,  1,  FWD,  1), true);
	}
}
