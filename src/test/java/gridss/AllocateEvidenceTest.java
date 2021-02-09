package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.sam.SAMFileUtil;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.FileHelper;
import com.google.common.collect.*;
import com.google.common.util.concurrent.MoreExecutors;
import gridss.cmdline.programgroups.Assembly;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.variant.variantcontext.FastGenotype;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class AllocateEvidenceTest extends IntermediateFilesTest {
	private void assertSymmetricalCalls(List<VariantContextDirectedEvidence> calls) {
		UnmodifiableIterator<VariantContextDirectedBreakpoint> it = Iterators.filter(calls.iterator(), VariantContextDirectedBreakpoint.class);
		assertSymmetrical(Lists.newArrayList(it));
	}
	private void assertSymmetrical(List<VariantContextDirectedBreakpoint> breakpoints) {
		Set<String> evidenceId = breakpoints.stream().map(bp -> bp.getEvidenceID()).collect(Collectors.toSet());
		Set<String> remoteEvidenceId = breakpoints.stream().map(bp -> bp.getRemoteEvidenceID()).collect(Collectors.toSet());
		assertEquals(evidenceId, remoteEvidenceId);
	}
	@Test
	public void should_annotate_reads() throws IOException {
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().minSize = 0;
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().minReads = 0;
		createInput(
				RP(0, 1, 10),
				DP(0, 1, "5M5S", true, 1, 10, "5M", true),
				DP(0, 2, "5M5S", true, 1, 10, "5M", true));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		FileHelper.copy(ses.getFile(), ses.getSVFile(), true);
		File assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		aes.ensureExtracted();
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(ses), ImmutableList.of(aes));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		List<VariantContextDirectedBreakpoint> vcfs = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedBreakpoint.class));
		assertEquals(2 * 1, vcfs.size()); // both breakends
		assertSymmetrical(vcfs);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(vcfs.iterator()), MoreExecutors.newDirectExecutorService()));
		assertSymmetricalCalls(results);
		assertEquals(vcfs.size(), results.size());
		VariantContextDirectedBreakpoint e = (VariantContextDirectedBreakpoint)results.get(0);
		assertEquals(2, e.getBreakpointEvidenceCount());
		assertEquals(2, e.getBreakendEvidenceCountSoftClip());
	}
	@Test
	public void should_apply_filters() throws IOException {
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().minSize = 0;
		createInput(
				RP(0, 1, 10),
				DP(0, 1, "5M5S", true, 1, 10, "5M", true));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		FileHelper.copy(ses.getFile(), ses.getSVFile(), true);
		File assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		aes.ensureExtracted();
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(ses), ImmutableList.of(aes));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		List<VariantContextDirectedBreakpoint> vcfs = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedBreakpoint.class));
		assertEquals(2 * 1, vcfs.size()); // both breakends
		assertSymmetrical(vcfs);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(vcfs.iterator()), MoreExecutors.newDirectExecutorService()));
		assertSymmetricalCalls(results);
		// single read support
		assertEquals(0, results.size());
	}
	@Test
	public void should_uniquely_assign() throws IOException, InterruptedException, ExecutionException {
		final int fragSize = 4;
		final int testSize = 16;
		final List<SAMRecord> in = new ArrayList<SAMRecord>();
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().writeFiltered = true;
		pc.getVariantCallingParameters().minScore = 0;
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
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.ASSEMBLY = ImmutableList.of(new File(testFolder.getRoot(), "assembly.bam"));
		createBAM(cmd.ASSEMBLY.get(0), AES().getHeader()); // to shut up the command-line parsing
		StubAssemblyEvidenceSource aes = new StubAssemblyEvidenceSource(pc, cmd.ASSEMBLY.get(0));
		aes.fragSize = fragSize;
		Collections.sort(ses.evidence, DirectedEvidenceOrder.ByNatural);
		createInput(in);
		VariantCaller vc = new VariantCaller(pc, ImmutableList.<SAMEvidenceSource>of(ses), ImmutableList.of(aes));
		ExecutorService threadpool = Executors.newSingleThreadExecutor();
		vc.callBreakends(output, threadpool);
		threadpool.shutdown();

		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		cmd.doWork(null);
		
		//List<IdsvVariantContext> annotated = getVcf(new File(testFolder.getRoot(), "out.vcf"), null);
		List<IdsvVariantContext> rawcalls = getVcf(output, null);
		List<IdsvVariantContext> calls = getVcf(cmd.OUTPUT_VCF, null);		
		assertSymmetrical(rawcalls.stream().map(x -> (VariantContextDirectedBreakpoint)x).collect(Collectors.toList()));
		assertSymmetrical(calls.stream().map(x -> (VariantContextDirectedBreakpoint)x).collect(Collectors.toList()));
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
		
		// each piece of evidence should be assigned to a single breakpoint so totals should match
		int rpCalls = calls.stream().mapToInt(v -> ((VariantContextDirectedBreakpoint)v).getBreakpointEvidenceCount()).sum();
		assertEquals(ses.evidence.size(), rpCalls);
		assertEquals(expectedEvidence, annotatedEvidence, 20); // floating point truncation on VCF is severe!
	}
	@Test
	public void should_filter_if_insufficient_reads() throws IOException {
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().minSize = 0;
		pc.getVariantCallingParameters().minReads = 3;
		createInput(
				RP(0, 1, 10),
				DP(0, 1, "5M5S", true, 1, 10, "5M", true));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		FileHelper.copy(ses.getFile(), ses.getSVFile(), true);
		File assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		aes.ensureExtracted();
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(ses), ImmutableList.of(aes));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		List<VariantContextDirectedBreakpoint> vcfs = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedBreakpoint.class));
		List<VariantContextDirectedEvidence> results = Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(vcfs.iterator()), MoreExecutors.newDirectExecutorService()));
		assertEquals(2, vcfs.size());
		assertEquals(0, results.size());
	}
	@Test
	public void insufficient_reads_filter_should_count_fragments() throws IOException {
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().minSize = 0;
		pc.getVariantCallingParameters().minReads = 2;
		createInput(
				RP(0, 1, 10),
				withReadName("read", DP(0, 1, "5M5S", true, 1, 10, "5M", true)),
				withReadName("read", DP(0, 2, "5M5S", true, 1, 10, "5M", true)),
				withReadName("read", DP(0, 1, "4M6S", true, 1, 10, "5M", true)),
				withReadName("read", DP(0, 1, "6M4S", true, 1, 10, "5M", true)));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		FileHelper.copy(ses.getFile(), ses.getSVFile(), true);
		File assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		aes.ensureExtracted();
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(ses), ImmutableList.of(aes));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		List<VariantContextDirectedBreakpoint> vcfs = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedBreakpoint.class));
		List<VariantContextDirectedEvidence> results = Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(vcfs.iterator()), MoreExecutors.newDirectExecutorService()));
		assertEquals(2, vcfs.size());
		assertEquals(0, results.size());
	}
	@Test
	public void should_filter_if_insufficient_quality() throws IOException {
		final ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().minScore = 1000;
		pc.getVariantCallingParameters().minSize = 0;
		pc.getVariantCallingParameters().minReads = 0;
		createInput(
				RP(0, 1, 10),
				DP(0, 1, "5M5S", true, 1, 10, "5M", true));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		FileHelper.copy(ses.getFile(), ses.getSVFile(), true);
		File assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		aes.ensureExtracted();
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(ses), ImmutableList.of(aes));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		List<VariantContextDirectedBreakpoint> vcfs = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedBreakpoint.class));
		List<VariantContextDirectedEvidence> results = Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(vcfs.iterator()), MoreExecutors.newDirectExecutorService()));
		assertEquals(0, results.size());
	}
	private void assembleSingleCategory(List<String> categories, File f, File ass, SAMRecord... reads) throws IOException {
		createBAM(f, SortOrder.coordinate, reads);
		ProcessingContext pc = getCommandlineContext(categories);
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, f, null, 0);
		SAMFileUtil.merge(ImmutableList.of(f), new File(pc.getFileSystemContext().getIntermediateDirectory(f), f.getName() + ".sv.bam"));
		ses.ensureMetrics();
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), ass);
		aes.assembleBreakends(MoreExecutors.newDirectExecutorService());
		// just want to index but can't find the htsjdk API for it
		SAMFileUtil.merge(ImmutableList.of(ass), new File(pc.getFileSystemContext().getIntermediateDirectory(ass), ass.getName() + ".sv.bam"));
		aes.ensureMetrics();
	}
	@Test
	public void should_support_batched_assembly() throws IOException, ExecutionException, InterruptedException {
		File n = new File(testFolder.getRoot(), "n.bam");
		File t = new File(testFolder.getRoot(), "t.bam");
		File a1 = new File(testFolder.getRoot(), "a1.bam");
		File a2 = new File(testFolder.getRoot(), "a2.bam");

		assembleSingleCategory(ImmutableList.of("Normal"), n, a1,
				withName("n1", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 1, "41M58S")))[0],
				withName("n2", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 1, "41M59S")))[0]);

		assembleSingleCategory(ImmutableList.of("Tumour"), t, a2,
				withName("t1", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 100, "41M58S")))[0],
				withName("t2", withSequence("AATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACCATTACCCGATCAAAACATATCAGAAATGATTGACGTATCACAAGCCGGAT", Read(0, 100, "41M59S")))[0]);

		ProcessingContext pc = getCommandlineContext();
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().minSize = 0;
		pc.getVariantCallingParameters().minReads = 0;
		SAMEvidenceSource fses1 = new SAMEvidenceSource(pc, n, null, 0);
		SAMEvidenceSource fses2 = new SAMEvidenceSource(pc, t, null, 1);
		AssemblyEvidenceSource faes1 = new AssemblyEvidenceSource(pc, ImmutableList.of(fses1, fses2), a1);
		AssemblyEvidenceSource faes2 = new AssemblyEvidenceSource(pc, ImmutableList.of(fses1, fses2), a2);
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(fses1, fses2), ImmutableList.of(faes1, faes2));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
		AllocateEvidence cmd = new AllocateEvidence();
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(faes1, faes2));
		cmd.setSamEvidenceSources(ImmutableList.of(fses1, fses2));
		cmd.OUTPUT_VCF = new File(testFolder.getRoot(), "annotated.vcf");
		cmd.doWork(MoreExecutors.newDirectExecutorService());
		List<VariantContext> results = getRawVcf(cmd.OUTPUT_VCF);

		for (VariantContext e : results) {
			for (Genotype g : e.getGenotypes()) {
				Object bsc = g.getExtendedAttribute("BSC");
				Object bassr = g.getExtendedAttribute("BASSR");
				Assert.assertEquals(bsc, bassr);
			}
		}
	}
}
