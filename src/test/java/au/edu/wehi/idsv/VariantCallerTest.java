package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;


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
		VariantCaller vc = new VariantCaller(pc, output, ImmutableList.<SAMEvidenceSource>of(ses), aes, null);
		vc.callBreakends(null);
		vc.annotateBreakpoints(null);
		List<IdsvVariantContext> annotated = getVcf(new File(testFolder.getRoot(), "out.vcf"), null);
		List<IdsvVariantContext> calls = getVcf(new File(new File(testFolder.getRoot(), "out.vcf.idsv.working"), "out.vcf.idsv.breakpoint.vcf"), null);
		// with no filtering, annotation should not change call set
		assertEquals(calls.size(), annotated.size());
		double expectedEvidence = 0;
		for (DirectedEvidence e : ses.evidence) {
			DiscordantReadPair bp = (DiscordantReadPair) e;
			expectedEvidence += bp.getBreakpointQual();
		}
		double annotatedEvidence = 0;
		for (IdsvVariantContext e : annotated) {
			annotatedEvidence += e.getPhredScaledQual();
		}
		// each piece of evidence should be assigned to a single breakpoint so totals should match
		assertEquals(expectedEvidence, annotatedEvidence, 20); // floating point truncation on VCF is severe!
	}
}
