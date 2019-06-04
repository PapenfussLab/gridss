package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;

import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import htsjdk.samtools.SAMRecord;

public class SequentialEvidenceAllocatorTest extends TestHelper {
	@Test
	public void should_uniquely_assign() {
		final int fragSize = 4;
		final int testSize = 64;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().breakendMargin = 0;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, null, 0, 0, fragSize);
		for (int i = 1; i < testSize; i++) {
			for (int j = 1; j < testSize; j++) {
				SAMRecord[] dp = withReadName(String.format("read-%d-%d", i, j), DP(0, i, "1M", true, 1, j, "1M", false));
				ses.evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
				ses.evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
				in.add(dp[0]);
				in.add(dp[1]);
			}
		}
		ses.evidence.sort(DirectedEvidenceOrder.ByNatural);
		AggregateEvidenceSource es = new AggregateEvidenceSource(pc, ImmutableList.of(ses), null);
		VariantCallIterator processor = new VariantCallIterator(es);
		ArrayList<VariantContextDirectedEvidence> calls = Lists.newArrayList(processor);
		calls.sort(VariantContextDirectedEvidence.ByBreakendStartEnd);
		
		SequentialEvidenceAllocator allocator = new SequentialEvidenceAllocator(pc, calls.iterator(), ses.evidence.iterator(), ImmutableList.<DirectedEvidence>of().iterator(), 4, true);
		ArrayList<VariantEvidenceSupport> result = Lists.newArrayList(allocator);
		
		assertEquals(ses.evidence.size(), result.stream().mapToInt(ves -> ves.support.size()).sum());
		for (VariantEvidenceSupport ves : result) {
			for (DirectedEvidence evidence : ves.support) {
				// check the assignment is valid
				assertTrue(ves.variant.getBreakendSummary().overlaps(evidence.getBreakendSummary()));
				// check the assignment is to the best option
				assertEquals(ves.variant.getPhredScaledQual(), calls.stream()
						.filter(v -> v.getBreakendSummary().overlaps(evidence.getBreakendSummary()))
						.mapToDouble(v -> v.getPhredScaledQual())
						.max().getAsDouble(), 0);
			}
		}
	}
	@Test
	public void should_assign_breakend_evidence_to_both_matching_breakpoint_and_matching_breakend() throws IOException, InterruptedException, ExecutionException {
		ArrayList<VariantContextDirectedEvidence> calls = new ArrayList<>();
		VariantContextDirectedEvidence be = (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 10), "")
				.phredScore(10).make();
		VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, BWD, 10, 1, BWD, 20), "")
				.phredScore(1).make();		
		calls.add(be);
		calls.add(bp);
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(getContext(), null, 0, 0, 100);
		ses.evidence.add(SCE(BWD, Read(0, 10, "10S10M")));
		SequentialEvidenceAllocator allocator = new SequentialEvidenceAllocator(getContext(), calls.iterator(), ses.evidence.iterator(), ImmutableList.<DirectedEvidence>of().iterator(), 1000000, true);
		ArrayList<VariantEvidenceSupport> result = Lists.newArrayList(allocator);
		VariantEvidenceSupport vesbp = result.stream().filter(ves -> ves.variant == bp).findFirst().orElse(null);
		VariantEvidenceSupport vesbe = result.stream().filter(ves -> ves.variant == be).findFirst().orElse(null);
		assertEquals(0, vesbe.support.size());
		assertEquals(1, vesbp.support.size());
	}
	@Test
	public void should_assign_breakend_evidence_to_breakend_if_no_matching_breakpoint() throws IOException, InterruptedException, ExecutionException {
		ArrayList<VariantContextDirectedEvidence> calls = new ArrayList<>();
		VariantContextDirectedEvidence be = (VariantContextDirectedEvidence)minimalBreakend()
				.breakend(new BreakendSummary(0, BWD, 10), "")
				.phredScore(10).make();
		calls.add(be);
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(getContext(), null, 0, 0, 100);
		ses.evidence.add(SCE(BWD, Read(0, 10, "10S10M")));
		SequentialEvidenceAllocator allocator = new SequentialEvidenceAllocator(getContext(), calls.iterator(), ses.evidence.iterator(), ImmutableList.<DirectedEvidence>of().iterator(), 1000000, true);
		ArrayList<VariantEvidenceSupport> result = Lists.newArrayList(allocator);
		assertEquals(1, result.get(0).support.size());
	}
	@Test
	public void RemoteOverlap_localLookup_should_split_on_referenceIndex_and_direction() {
		final ProcessingContext pc = getContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().breakendMargin = 0;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, null, 0, 0, 100);
		ses.evidence.add(SCE(FWD, Read(0, 1, "1S1M1S")));
		ses.evidence.add(SCE(BWD, Read(0, 1, "1S1M1S")));
		//ses.evidence.add(SCE(BWD, Read(1, 1, "1S1M1S")));
		//ses.evidence.add(SCE(FWD, Read(1, 1, "1S1M1S")));
		SequentialEvidenceAllocator allocator = new SequentialEvidenceAllocator(pc,
				ImmutableList.of((VariantContextDirectedBreakpoint)TestHelper.minimalVariant()
						.breakpoint(new BreakpointSummary(1, FWD, 1, 1, 1, 0, FWD, 1, 1, 1), "")
						.phredScore(1)
						.id("call")
						.make()).iterator(),
				ses.evidence.iterator(),
				ImmutableList.<DirectedEvidence>of().iterator(),
				Integer.MAX_VALUE, true);
		// shouldn't crash
		VariantEvidenceSupport ves = allocator.next();
		Assert.assertTrue(ves.support.stream().allMatch(e -> ves.variant.getBreakendSummary().overlaps(e.getBreakendSummary())));
	}
}
