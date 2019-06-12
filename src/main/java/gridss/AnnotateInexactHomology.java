package gridss;

import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.alignment.BreakpointHomology;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;

import java.util.Iterator;
import java.util.concurrent.ExecutorService;

public class AnnotateInexactHomology extends VcfTransformCommandLineProgram {
	@Override
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		Iterator<VariantContextDirectedEvidence> it = new ParallelTransformIterator<VariantContextDirectedEvidence, VariantContextDirectedEvidence>(
				calls,
				call -> (call instanceof VariantContextDirectedBreakpoint) ? BreakpointHomology.annotate(getContext(), (VariantContextDirectedBreakpoint)call) : call,
				WORKER_THREADS + 1,
				threadpool);
		return new AutoClosingIterator<>(it, calls);
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateInexactHomology().instanceMain(argv));
    }
}
