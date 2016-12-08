package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class MaximalEvidenceCliqueIteratorTest extends TestHelper {
	@Test
	public void breakends_should_have_matching_event_id() {
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new MaximalEvidenceCliqueIterator(getContext(), ImmutableList.<DirectedEvidence>of(NRRP(DP(0, 1, "1M", true, 0, 500, "1M", false))).iterator(), FWD, BWD, new SequentialIdGenerator("test")));
		assertEquals(2, calls.size());
		assertEquals("test1", calls.get(0).getAttribute("EVENT"));
	}
	@Test
	public void breakends_mateids_should_match() {
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new MaximalEvidenceCliqueIterator(getContext(), ImmutableList.<DirectedEvidence>of(NRRP(DP(0, 1, "1M", true, 0, 500, "1M", false))).iterator(), FWD, BWD, new SequentialIdGenerator("test")));
		assertEquals(2, calls.size());
		assertEquals(calls.get(0).getID(), calls.get(1).getAttribute("PARID"));
		assertEquals(calls.get(1).getID(), calls.get(0).getAttribute("PARID"));
	}
	@Test
	public void margin_should_expand_and_contract_past_chromosome_end() {
		BreakpointSummary circularbp = new BreakpointSummary(0, BWD, 1, 0, FWD, POLY_A.length);
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new MaximalEvidenceCliqueIterator(getContext(), ImmutableList.<DirectedEvidence>of(new MockDirectedBreakpoint(circularbp)).iterator(), BWD, FWD, new SequentialIdGenerator("test")));
		assertEquals(2, calls.size());
		assertEquals(circularbp, calls.get(0).getBreakendSummary());
		assertEquals(circularbp.remoteBreakpoint(), calls.get(1).getBreakendSummary());
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_require_valid_evidence() {
		BreakpointSummary outOfBounds = new BreakpointSummary(0, BWD, -10, 0, FWD, 12000);
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new MaximalEvidenceCliqueIterator(getContext(), ImmutableList.<DirectedEvidence>of(new MockDirectedBreakpoint(outOfBounds)).iterator(), BWD, FWD, new SequentialIdGenerator("test")));
		assertEquals(2, calls.size());
	}
}
