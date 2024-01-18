package au.edu.wehi.idsv;

import au.edu.wehi.idsv.util.ErrorIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.QueryInterval;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class VariantCallIteratorTest extends IntermediateFilesTest {
	@Test
	public void margin_should_expand_and_contract_past_chromosome_end()  {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 1, 0, FWD, POLY_A.length)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertTrue(result.get(0) instanceof VariantContextDirectedBreakpoint);
		assertEquals(new BreakpointSummary(0, BWD, 1, 0, FWD, POLY_A.length), result.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(0, FWD, POLY_A.length, 0, BWD, 1), result.get(1).getBreakendSummary());
	}
	@Test
	public void evidence_should_not_expand()  {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(SR(withSequence("TTTT", Read(0, 10, "1M3S"))[0], withSequence("TTT", Read(1, 10, "3M"))[0]));
		list.add(SR(withSequence("TTTT", Read(0, 11, "1M3S"))[0], withSequence("TTT", Read(1, 11, "3M"))[0]));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4, result.size());
	}
	@Test
	public void should_call_maximal_clique()  throws InterruptedException {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 15, 25, 1, BWD, 35, 35, 45)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertEquals(new BreakpointSummary(0, FWD, 17, 15, 20, 1, BWD, 37, 35, 40), result.get(0).getBreakendSummary());
	}
	@Test
	public void singleton_should_call_both_breakends() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 15, 10, 20, 1, BWD, 30, 30, 40)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
		assertTrue(result.get(0) instanceof VariantContextDirectedBreakpoint);
		assertEquals(new BreakpointSummary(0, FWD, 15, 10, 20, 1, BWD, 35, 30, 40), result.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(1, BWD, 35, 30, 40, 0, FWD, 15, 10, 20), result.get(1).getBreakendSummary());
	}
	@Test(expected=RuntimeException.class)
	public void should_rethrow_worker_thread_exceptions() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, -1, -1, -1, 1, BWD, 30, 30, 40)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		Lists.newArrayList(ecp);
	}
	@Test
	public void should_call_orientations()  {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, FWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, FWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 10, 20, 1, FWD, 30, 30, 40)));
		list.add(new MockDirectedBreakpoint(new BreakpointSummary(0, BWD, 10, 10, 20, 1, BWD, 30, 30, 40)));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(4 * 2, result.size());
	}
	@Test
	public void interval_caller_should_filter_calls_in_which_neither_breakend_starts_in_interval()  {
		createInput(
				RP(0, 1, 2, 1),
				DP(0, 1, "1M", true, 1, 100, "1M", false));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		VariantCallIterator ecp = new VariantCallIterator(
				new AggregateEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), null, SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition),
				new QueryInterval[] {new QueryInterval(0, 1, 10) },
				0);
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
	}
	@Test
	public void should_call_breakends() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		list.add(SCE(FWD, Read(0, 1, "10M10S")));
		list.add(SCE(BWD, Read(0, 100, "10S10M")));
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(2, result.size());
	}
	@Test
	public void add_to_deque_should_block() {
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		for (int i = 0; i < 1024; i++) {
			list.add(SCE(FWD, Read(0, i, "10M10S")));
		}
		VariantCallIterator ecp = new VariantCallIterator(getContext(), list.iterator());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
		assertEquals(1024, result.size());
	}

	/**
	 * https://github.com/PapenfussLab/gridss/issues/267
	 */
	@Test(expected = RuntimeException.class)
	public void iterator_exception_should_surface_275() {
		VariantCallIterator ecp = new VariantCallIterator(getContext(), new ErrorIterator<>());
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(ecp);
	}
	/**
	 * https://github.com/PapenfussLab/gridss/issues/267
	 */
	@Test(expected = RuntimeException.class)
	public void iterator_exception_should_surface_from_parsing_error_275() throws IOException {
		ArrayList<SAMEvidenceSource> sess = Lists.newArrayList(new MockSAMEvidenceSource(getContext(), new File("src/test/resources/malformedsv.sam")));
		MockAssemblyEvidenceSource aes = new MockAssemblyEvidenceSource(getContext(), sess, new File("src/test/resources/empty.bam"));
		VariantCaller vc = new VariantCaller(getContext(), sess, ImmutableList.of(aes));
		ExecutorService threadpool = Executors.newFixedThreadPool(2, new ThreadFactoryBuilder().setDaemon(true).setNameFormat("Test275-%d").build());
		vc.callBreakends(new File(testFolder.getRoot(), "out275.vcf"), threadpool);
	}
}
