package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import htsjdk.samtools.QueryInterval;

public class VariantCallIteratorTest extends IntermediateFilesTest {
	@Test
	public void margin_should_expand_and_contract_past_chromosome_end() throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 1, 0, FWD, POLY_A.length)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertTrue(result.get(0) instanceof VariantContextDirectedBreakpoint);
		assertEquals(new BreakpointSummary(0, BWD, 1, 0, FWD, POLY_A.length), result.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(0, FWD, POLY_A.length, 0, BWD, 1), result.get(1).getBreakendSummary());
	}
	@Test
	public void evidence_should_not_expand() throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(SR(withSequence("TTTT", Read(0, 10, "1M3S"))[0], withSequence("TTT", Read(1, 10, "3M"))[0]));
		list.add(SR(withSequence("TTTT", Read(0, 11, "1M3S"))[0], withSequence("TTT", Read(1, 11, "3M"))[0]));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4, result.size());
	}
	@Test
	public void should_call_maximal_clique()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 15, 25, 1, BWD, 35, 35, 45)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertEquals(new BreakpointSummary(0, FWD, 17, 15, 20, 1, BWD, 37, 35, 40), result.get(0).getBreakendSummary());
	}
	@Test
	public void singleton_should_call_both_breakends()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 10, 20, 1, BWD, 30, 30, 40)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertTrue(result.get(0) instanceof VariantContextDirectedBreakpoint);
		assertEquals(new BreakpointSummary(0, FWD, 15, 10, 20, 1, BWD, 35, 30, 40), result.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(1, BWD, 35, 30, 40, 0, FWD, 15, 10, 20), result.get(1).getBreakendSummary());
	}
	@Test(expected=RuntimeException.class)
	public void should_rethrow_worker_thread_exceptions()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, -1, -1, -1, 1, BWD, 30, 30, 40)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		Lists.newArrayList(ecp);
	}
	@Test
	public void should_call_orientations()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, FWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 10, 20, 1, FWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4 * 2, result.size());
	}
	@Test
	public void interval_caller_should_filter_calls_in_which_neither_breakend_starts_in_interval()  throws InterruptedException {
		createInput(
				RP(0, 1, 2, 1),
				DP(0, 1, "1M", true, 1, 100, "1M", false));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		VariantCallIterator ecp = new VariantCallIterator(
				new AggregateEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), null),
				new QueryInterval[] {new QueryInterval(0, 1, 10) },
				0);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
	}
	@Test
	public void should_call_breakends() throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(SCE(FWD, Read(0, 1, "10M10S")));
		list.add(SCE(FWD, Read(0, 2, "9M10S")));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
	}
}
