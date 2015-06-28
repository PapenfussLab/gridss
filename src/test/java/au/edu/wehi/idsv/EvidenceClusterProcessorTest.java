package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

public class EvidenceClusterProcessorTest extends TestHelper {
	@Test
	public void margin_should_expand_and_contract_past_chromosome_end() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 1, 1, 0, FWD, POLY_A.length, POLY_A.length)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertTrue(result.get(0) instanceof VariantContextDirectedBreakpoint);
		assertEquals(new BreakpointSummary(0, BWD, 1, 1, 0, FWD, POLY_A.length, POLY_A.length), result.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(0, FWD, POLY_A.length, POLY_A.length, 0, BWD, 1, 1), result.get(1).getBreakendSummary());
	}
	@Test
	public void evidence_should_expand_by_margin_then_shrink() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(SCE(FWD, withSequence("TTTT", Read(0, 10, "1M3S"))[0], withSequence("TTTT", Read(1, 10, "4M"))[0]));
		list.add(SCE(FWD, withSequence("TTTT", Read(0, 13, "1M3S"))[0], withSequence("TTTT", Read(1, 13, "4M"))[0]));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		BreakpointSummary bp = ((VariantContextDirectedBreakpoint)result.get(0)).getBreakendSummary();
		// expand by 3bp, overlap is 10-13
		// then shrink back down to call the middle of the interval
		assertEquals(11, bp.start);
		assertEquals(11, bp.end);
	}
	@Test
	public void should_call_maximal_clique() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 20, 1, BWD, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 25, 1, BWD, 35, 45)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertEquals(new BreakpointSummary(0, FWD, 15, 20, 1, BWD, 35, 40), result.get(0).getBreakendSummary());
	}
	@Test
	public void singleton_should_call_both_breakends() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 20, 1, BWD, 30, 40)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertTrue(result.get(0) instanceof VariantContextDirectedBreakpoint);
		assertEquals(new BreakpointSummary(0, FWD, 10, 20, 1, BWD, 30, 40), result.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(1, BWD, 30, 40, 0, FWD, 10, 20), result.get(1).getBreakendSummary());
	}
	@Test(expected=RuntimeException.class)
	public void should_rethrow_worker_thread_exceptions() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, -1, -1, 1, BWD, 30, 40)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(getContext(), list.iterator());
		Lists.newArrayList(ecp);
	}
	@Test
	public void should_call_orientations() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 20, 1, FWD, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 20, 1, BWD, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 20, 1, FWD, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 20, 1, BWD, 30, 40)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4 * 2, result.size());
	}
}
