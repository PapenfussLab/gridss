package gridss;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ExecutionException;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.VariantCaller;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;

public class AllocateEvidenceTest extends IntermediateFilesTest {
	@Test
	public void should_annotate_reads() throws IOException {
		createInput(
				RP(0, 1, 10),
				DP(0, 1, "5M5S", true, 1, 10, "4S5M1D5M6S", true));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, 0);
		ses.ensureMetrics();
		File assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getContext(), ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		VariantCaller caller = new VariantCaller(getContext(), ImmutableList.of(ses), aes);
		caller.callBreakends(output, null);
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(getContext());
		cmd.setAssemblySource(aes);
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		List<VariantContextDirectedBreakpoint> vcfs = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedBreakpoint.class));
		List<VariantContextDirectedBreakpoint> results = Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(vcfs.iterator()), MoreExecutors.newDirectExecutorService()));
		assertEquals(vcfs.size(), results.size());
	}
	@Test
	public void should_uniquely_assign() throws IOException, InterruptedException, ExecutionException {
		final int fragSize = 4;
		final int testSize = 64;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		StubSAMEvidenceSource ses = new StubSAMEvidenceSource(pc, input, 0, 0, fragSize);
		for (int i = 1; i < testSize; i++) {
			for (int j = 1; j < testSize; j++) {
				SAMRecord[] dp = withReadName(String.format("read-%d-%d", i, j), DP(0, i, "1M", true, 1, j, "1M", false));
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
		vc.callBreakends(output, MoreExecutors.newDirectExecutorService());
		
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(aes);
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		cmd.ASSEMBLY = new File(testFolder.getRoot(), "assembly.bam");
		createBAM(cmd.ASSEMBLY, SortOrder.coordinate, Collections.emptyList()); // to shut up the command-line parsing
		cmd.doWork(null);
		
		//List<IdsvVariantContext> annotated = getVcf(new File(testFolder.getRoot(), "out.vcf"), null);
		List<IdsvVariantContext> rawcalls = getVcf(output, null);
		List<IdsvVariantContext> calls = getVcf(cmd.OUTPUT_VCF, null);
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
		double rawEvidence = rawcalls.stream().mapToDouble(e -> e.getPhredScaledQual()).sum();
		// unique assignment must result in the evidence total being, at most, the same 
		assertTrue(annotatedEvidence <= rawEvidence);
		
		sanityCheck(rawcalls, calls);
		// each piece of evidence should be assigned to a single breakpoint so totals should match
		assertEquals(expectedEvidence, annotatedEvidence, 20); // floating point truncation on VCF is severe!
	}
	@Test
	public void should_call_cliques() throws IOException, InterruptedException, ExecutionException {
		should_call_cliques(0);
	}
	@Test
	public void reallocation_margin_should_match_max_clique_calculation_margin() throws IOException, InterruptedException, ExecutionException {
		should_call_cliques(1);
	}
	private void should_call_cliques(int breakendMargin) throws IOException, InterruptedException, ExecutionException {
		final int fragSize = 4;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		pc.getVariantCallingParameters().minScore = 0.001;
		pc.getVariantCallingParameters().breakendMargin = breakendMargin;
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
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(aes);
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		cmd.ASSEMBLY = new File(testFolder.getRoot(), "assembly.bam");
		createBAM(cmd.ASSEMBLY, SortOrder.coordinate, Collections.emptyList()); // to shut up the command-line parsing
		cmd.doWork(null);
		
		List<IdsvVariantContext> rawcalls = getVcf(output, null);
		List<IdsvVariantContext> calls = getVcf(cmd.OUTPUT_VCF, null);
		
		sanityCheck(rawcalls, calls);
		assertEquals(3 * ((DirectedBreakpoint)ses.evidence.get(0)).getBreakpointQual(), calls.get(0).getPhredScaledQual(), 0.01);
	}
	private void sanityCheck(List<IdsvVariantContext> rawcalls, List<IdsvVariantContext> calls) {
		for (IdsvVariantContext raw : rawcalls) {
			Optional<IdsvVariantContext> ann = calls.stream().filter(v -> v.getID().equals(raw.getID())).findFirst();
			if (ann.isPresent()) {
				// support can never go up after reallocation
				assertTrue(raw.getPhredScaledQual() >= ann.get().getPhredScaledQual());
			}
		}
	}
}
