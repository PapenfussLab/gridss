package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import htsjdk.samtools.SAMRecord;

public class SequentialEvidenceAllocatorTest extends TestHelper {
	@Test
	public void should_uniquely_assign() throws IOException, InterruptedException, ExecutionException {
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
		ExecutorService threadpool = Executors.newFixedThreadPool(1);
		EvidenceClusterProcessor processor = new EvidenceClusterProcessor(threadpool, es);
		ArrayList<VariantContextDirectedEvidence> calls = Lists.newArrayList(processor);
		threadpool.shutdown();
		calls.sort(VariantContextDirectedEvidence.ByBreakendStartEnd);
		
		SequentialEvidenceAllocator allocator = new SequentialEvidenceAllocator(pc, calls.iterator(), ses.evidence.iterator(), 4, true);
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
}
