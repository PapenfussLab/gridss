package gridss;

import java.util.concurrent.ExecutorService;

import au.edu.wehi.idsv.SequentialCoverageAnnotator;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;

public class AnnotateReferenceCoverage extends VcfTransformCommandLineProgram {
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls, ExecutorService threadpool) {
		return new SequentialCoverageAnnotator<VariantContextDirectedBreakpoint>(getContext(), getSamEvidenceSources(), calls, threadpool);
	}
}
