package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.util.concurrent.MoreExecutors;

import htsjdk.samtools.SAMRecord;


public class VariantCallerTest extends IntermediateFilesTest {
	@Test
	public void shouldCallMaxCliques() {
		final int fragSize = 4;
		final int testSize = 64;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, input, 0, 0, fragSize);
		for (int i = 1; i < testSize; i++) {
			for (int j = 1; j < testSize; j++) {
				SAMRecord[] dp = DP(0, i, "1M", true, 1, j, "1M", false);
				ses.evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
				ses.evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
				in.add(dp[0]);
				in.add(dp[1]);
			}
		}
		StubAssemblyEvidenceSource aes = new StubAssemblyEvidenceSource(pc);
		aes.fragSize = fragSize;
		Collections.sort(ses.evidence, DirectedEvidenceOrder.ByNatural);
		createInput(in);
		VariantCaller vc = new VariantCaller(pc, ImmutableList.<SAMEvidenceSource>of(ses), aes);
		ExecutorService threadpool = Executors.newSingleThreadExecutor();
		vc.callBreakends(output, threadpool);
		threadpool.shutdown();
		//List<IdsvVariantContext> annotated = getVcf(new File(testFolder.getRoot(), "out.vcf"), null);
		List<IdsvVariantContext> calls = getVcf(output, null);
		// with no filtering, annotation should not change call set
		double expectedEvidence = 0;
		for (DirectedEvidence e : ses.evidence) {
			DiscordantReadPair bp = (DiscordantReadPair) e;
			expectedEvidence += bp.getBreakpointQual();
		}
		double annotatedEvidence = 0;
		for (IdsvVariantContext e : calls) {
			annotatedEvidence += e.getPhredScaledQual();
		}
		// each piece of evidence should be assigned to a single breakpoint so totals should match
		// TODO: what's the total for these overlapping cliques?
		//assertEquals(expectedEvidence, annotatedEvidence, 20); // floating point truncation on VCF is severe!
	}
	@Test
	public void should_call_cliques() {
		final int fragSize = 4;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		pc.getVariantCallingParameters().breakendMargin = 0;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, input, 0, 0, fragSize);
		for (int i = 1; i <= 5; i++) {
			SAMRecord[] dp = DP(0, i, "1M", true, 1, i, "1M", true);
			ses.evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
			ses.evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
			in.add(dp[0]);
			in.add(dp[1]);
		}
		StubAssemblyEvidenceSource aes = new StubAssemblyEvidenceSource(pc);
		aes.fragSize = fragSize;
		Collections.sort(ses.evidence, DirectedEvidenceOrder.ByNatural);
		createInput(in);
		VariantCaller vc = new VariantCaller(pc, ImmutableList.<SAMEvidenceSource>of(ses), aes);
		vc.callBreakends(output, MoreExecutors.newDirectExecutorService());
		//vc.annotateBreakpoints(MoreExecutors.newDirectExecutorService());
		//List<IdsvVariantContext> annotated = getVcf(new File(testFolder.getRoot(), "out.vcf"), null);
		List<IdsvVariantContext> calls = getVcf(output, null);
		// start/end ramps are not max cliques but the rest are
		assertEquals(2 * 3, calls.size());
		for (IdsvVariantContext variant : calls) {
			assertEquals(3 * ((DirectedBreakpoint)ses.evidence.get(0)).getBreakpointQual(), variant.getPhredScaledQual(), 0.01);
		}
	}
}
