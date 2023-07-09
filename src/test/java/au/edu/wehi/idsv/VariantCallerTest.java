package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import com.google.common.util.concurrent.MoreExecutors;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static org.junit.Assert.assertEquals;


public class VariantCallerTest extends IntermediateFilesTest {
	@Test
	public void shouldCallMaxCliques() throws IOException {
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
		VariantCaller vc = new VariantCaller(pc, ImmutableList.<SAMEvidenceSource>of(ses), ImmutableList.of(aes));
		ExecutorService threadpool = Executors.newSingleThreadExecutor();
		vc.callBreakends(output, threadpool, 0, 1);
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
	public void should_call_cliques() throws IOException {
		final int fragSize = 4;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		pc.getVariantCallingParameters().breakendMargin = 0;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, input, 0, 0, fragSize);
		List<IdsvVariantContext> calls = go_call_simple_cliques(fragSize, in, pc, ses);
		// start/end ramps are not max cliques but the rest are
		assertEquals(2 * 3, calls.size());
		for (IdsvVariantContext variant : calls) {
			assertEquals(3 * ((DirectedBreakpoint)ses.evidence.get(0)).getBreakpointQual(), variant.getPhredScaledQual(), 0.01);
		}
	}

	private List<IdsvVariantContext> go_call_simple_cliques(int fragSize, List<SAMRecord> in, ProcessingContext pc, StubSAMEvidenceSource ses) throws IOException {
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
		VariantCaller vc = new VariantCaller(pc, ImmutableList.<SAMEvidenceSource>of(ses), ImmutableList.of(aes));
		vc.callBreakends(output, MoreExecutors.newDirectExecutorService(), 0,1);
		//vc.annotateBreakpoints(MoreExecutors.newDirectExecutorService());
		//List<IdsvVariantContext> annotated = getVcf(new File(testFolder.getRoot(), "out.vcf"), null);
		return getVcf(output, null);
	}

	@Test
	public void should_drop_low_qual() throws IOException {
		final int fragSize = 4;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = false;
		pc.getVariantCallingParameters().minScore = 1000;
		pc.getVariantCallingParameters().breakendMargin = 0;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, input, 0, 0, fragSize);
		List<IdsvVariantContext> calls = go_call_simple_cliques(fragSize, in, pc, ses);
		// same as should_call_cliques test
		assertEquals(0, calls.size());
	}

	/**
	 * These are likely to be false positives caused by the reads slightly outside of the
	 * 99% fragment size distribution
	 */
	@Test
	public void should_drop_imprecise_small_deletions() throws IOException {
		for (boolean writeFiltered : new boolean[] { false, true}) {
			ProcessingContext pc = getCommandlineContext();
			pc.getVariantCallingParameters().writeFiltered = writeFiltered;
			pc.getVariantCallingParameters().minScore = 1;
			int fragSize = 500;
			StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, input, 0, 0, 500);
			for (int i = fragSize + 1; i < fragSize + 200; i++) {
				SAMRecord[] dp = DP(0, 1, "100M", true, 0, i + 100, "100M", false);
				ses.evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
				ses.evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
			}
			// Single breakend
			ses.evidence.add(SCE(FWD, ses, Read(1, 1, "50M50S")));
			ses.evidence.add(SCE(FWD, ses, Read(1, 2, "49M50S")));
			ses.evidence.add(SCE(FWD, ses, Read(1, 3, "48M50S")));
			// exact small deletion
			ses.evidence.add(SR(ses, Read(2, 110, "50S50M"), Read(2, 200, "50M")));
			SAMRecord[] dp = DP(2, 100, "100M", true, 2, 100+500+200-100, "100M", false);
			ses.evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
			ses.evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
			Collections.sort(ses.evidence, DirectedEvidenceOrder.ByNatural);
			VariantCaller vc = new VariantCaller(pc, ImmutableList.<SAMEvidenceSource>of(ses), ImmutableList.of());
			vc.callBreakends(output, MoreExecutors.newDirectExecutorService(), 0, 1);
			List<IdsvVariantContext> result = getVcf(output, null);
			if (writeFiltered) {
				Assert.assertTrue(result.size() > 3);
			} else {
				Assert.assertEquals(3, result.size());
				Assert.assertTrue(result.stream().allMatch(v -> !v.hasAttribute("IMPRECISE")));
			}
			output.delete();
		}
	}
}
