
package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;

public class BreakendMaximalEvidenceCliqueIteratorTest extends TestHelper {
	@Test
	public void should_ignore_breakpoint_evidence() {
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new BreakendMaximalEvidenceCliqueIterator(
				getContext(),
				ImmutableList.<DirectedEvidence>of(NRRP(DP(0, 1, "1M", true, 0, 500, "1M", false))).iterator(),
				FWD,
				new SequentialIdGenerator("test")));
		assertEquals(0, calls.size());
	}
	@Test
	public void should_call_breakend_evidence() {
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new BreakendMaximalEvidenceCliqueIterator(
				getContext(),
				ImmutableList.<DirectedEvidence>of(SCE(FWD, Read(0, 1, "50M50S"))).iterator(),
				FWD,
				new SequentialIdGenerator("test")));
		assertEquals(1, calls.size());
	}
	@Test
	public void should_call_correct_direction() {
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new BreakendMaximalEvidenceCliqueIterator(
				getContext(),
				ImmutableList.<DirectedEvidence>of(SCE(FWD, Read(0, 1, "50M50S"))).iterator(),
				BWD,
				new SequentialIdGenerator("test")));
		assertEquals(0, calls.size());
	}
	@Test
	public void should_call_score_all_variants() {
		ImmutableList<DirectedEvidence> evidence = ImmutableList.<DirectedEvidence>of(
				SCE(FWD, Read(0, 1, "50M50S")),
				SCE(FWD, Read(0, 1, "50M49S")),
				SCE(FWD, Read(0, 100, "50M49S"))
				);
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new BreakendMaximalEvidenceCliqueIterator(
				getContext(),
				evidence.iterator(),
				FWD,
				new SequentialIdGenerator("test")));
		assertEquals(2, calls.size());
		assertEquals(evidence.get(0).getBreakendQual() + evidence.get(1).getBreakendQual(), calls.get(0).getPhredScaledQual(), 1e-6);
		assertEquals(evidence.get(2).getBreakendQual(), calls.get(1).getPhredScaledQual(), 1e-6);
	}
	@Test
	public void should_call_overlapping_intervals() {
		ImmutableList<DirectedEvidence> evidence = ImmutableList.<DirectedEvidence>of(
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 5, 5, 10)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 7, 7, 9)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 10, 10, 10)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 11, 11, 11)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 100, 100, 200)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 110, 110, 190)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 120, 120, 180)),
				new MockDirectedEvidence(new BreakendSummary(0, FWD, 200, 200, 200))
				);
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new BreakendMaximalEvidenceCliqueIterator(
				getContext(),
				evidence.iterator(),
				FWD,
				new SequentialIdGenerator("test")));
		assertEquals(new BreakendSummary(0, FWD, 8, 7, 9), calls.get(0).getBreakendSummary());
		assertEquals(new BreakendSummary(0, FWD, 10, 10, 10), calls.get(1).getBreakendSummary());
		assertEquals(new BreakendSummary(0, FWD, 11, 11, 11), calls.get(2).getBreakendSummary());
		assertEquals(evidence.get(3).getBreakendQual(), calls.get(2).getPhredScaledQual(), 1e-6);
		assertEquals(new BreakendSummary(0, FWD, 150, 120, 180), calls.get(3).getBreakendSummary());
		assertEquals(new BreakendSummary(0, FWD, 200, 200, 200), calls.get(4).getBreakendSummary());
	}
}
