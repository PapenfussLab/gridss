package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;

public class EvidenceClusterProcessorTest extends TestHelper {
	@Test
	public void margin_should_expand_and_contract_past_chromosome_end() throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 1, 0, FWD, POLY_A.length)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(MoreExecutors.directExecutor(), getContext(), list);
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
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(MoreExecutors.directExecutor(), getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4, result.size());
	}
	@Test
	public void should_call_maximal_clique()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 15, 25, 1, BWD, 35, 35, 45)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(MoreExecutors.directExecutor(), getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertEquals(new BreakpointSummary(0, FWD, 17, 15, 20, 1, BWD, 37, 35, 40), result.get(0).getBreakendSummary());
	}
	@Test
	public void singleton_should_call_both_breakends()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 10, 20, 1, BWD, 30, 30, 40)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(MoreExecutors.directExecutor(), getContext(), list);
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
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(MoreExecutors.directExecutor(), getContext(), list);
		Lists.newArrayList(ecp);
	}
	@Test
	public void should_call_orientations()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, FWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 10, 20, 1, FWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(MoreExecutors.directExecutor(), getContext(), list);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4 * 2, result.size());
	}
	@Test
	public void should_run_in_parallel() throws InterruptedException {
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(getContext(), null, 0, 100, 300);
		for (int i = 1; i < POLY_A.length; i++) {
			for (int j = 0; j < 2; j++) {
				ses.evidence.add(NRRP(ses, DP(j, i, "1M", true, j, i, "1M", true)));
				ses.evidence.add(NRRP(ses, DP(j, i, "1M", true, j, i, "1M", false)));
				ses.evidence.add(NRRP(ses, DP(j, i, "1M", false, j, i, "1M", true)));
				ses.evidence.add(NRRP(ses, DP(j, i, "1M", false, j, i, "1M", false)));
			}
		}
		ExecutorService threadpool = Executors.newFixedThreadPool(8);
		EvidenceClusterProcessor ecp = new EvidenceClusterProcessor(threadpool, new AggregateEvidenceSource(getContext(), ImmutableList.of(ses), null));
		threadpool.shutdown();
		// tasks should be blocking on the output queue because we haven't pulled anything from it yet
		assertFalse(threadpool.awaitTermination(10, TimeUnit.MILLISECONDS));
		Lists.newArrayList(ecp); // traverse
		// tasks should now all be complete
		assertTrue(threadpool.awaitTermination(10, TimeUnit.MILLISECONDS));
	}
}
