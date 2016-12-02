package gridss;

import java.util.Iterator;
import java.util.concurrent.ExecutorService;

import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.alignment.BreakpointHomology;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;

public class AnnotateInexactHomology extends VcfTransformCommandLineProgram {
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls, ExecutorService threadpool) {
		Iterator<VariantContextDirectedBreakpoint> it = new ParallelTransformIterator<VariantContextDirectedBreakpoint, VariantContextDirectedBreakpoint>(
				calls, call -> BreakpointHomology.annotate(getContext(), call), WORKER_THREADS + 1, threadpool);
		return new AutoClosingIterator<>(it, calls);
	}
}
