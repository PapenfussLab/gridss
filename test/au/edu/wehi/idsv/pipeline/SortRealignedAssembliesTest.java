package au.edu.wehi.idsv.pipeline;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;

import java.io.IOException;
import java.util.EnumSet;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;

import com.google.common.collect.ImmutableList;


public class SortRealignedAssembliesTest extends IntermediateFilesTest {
	public VariantContextDirectedBreakpoint BP(String id, BreakpointSummary bp) {
		IdsvVariantContextBuilder builder = minimalBreakend()
				.breakpoint(bp, null);
		builder.id(id);
		return (VariantContextDirectedBreakpoint)builder.make();
	}
	@Test
	public void should_sort_by_remote_breakend_coordinates() throws IOException {
		ProcessingContext pc = getCommandlineContext(false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.<SAMEvidenceSource>of(), output);
		createVCF(pc, pc.getFileSystemContext().getAssemblyVcf(output),
				BP("1", new BreakpointSummary(0, FWD, 1, 1, 2, FWD, 10, 10)),
				BP("2", new BreakpointSummary(1, FWD, 1, 1, 1, FWD, 10, 10)),
				BP("3", new BreakpointSummary(1, FWD, 5, 5, 1, FWD, 1, 1)),
				BP("4", new BreakpointSummary(1, FWD, 10, 10, 0, FWD, 1, 1)),
				BP("5", new BreakpointSummary(2, FWD, 1, 1, 1, FWD, 2, 2))
				);
		pc.getFileSystemContext().getRealignmentFastq(output).createNewFile();
		createBAM(pc.getFileSystemContext().getRealignmentBam(output), SortOrder.unsorted);
		SortRealignedAssemblies sra = new SortRealignedAssemblies(pc, aes);
		sra.process(EnumSet.of(ProcessStep.SORT_REALIGNED_ASSEMBLIES));
		sra.close();
		
		List<IdsvVariantContext> result = getRA(aes);
		assertEquals(5, result.size());
		assertEquals("4", result.get(0).getID());
		assertEquals("3", result.get(1).getID());
		assertEquals("5", result.get(2).getID());
		assertEquals("2", result.get(3).getID());
		assertEquals("1", result.get(4).getID());
	}
	@Test
	public void should_ignore_breakend() throws IOException {
		ProcessingContext pc = getCommandlineContext(false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.<SAMEvidenceSource>of(), output);
		createVCF(pc, pc.getFileSystemContext().getAssemblyVcf(output),
				minimalBreakend().make()
				);
		pc.getFileSystemContext().getRealignmentFastq(output).createNewFile();
		createBAM(pc.getFileSystemContext().getRealignmentBam(output), SortOrder.unsorted); // mock realignment
		SortRealignedAssemblies sra = new SortRealignedAssemblies(pc, aes);
		sra.process(EnumSet.of(ProcessStep.SORT_REALIGNED_ASSEMBLIES));
		sra.close();
		
		List<IdsvVariantContext> result = getRA(aes);
		assertEquals(0, result.size());
	}
}
