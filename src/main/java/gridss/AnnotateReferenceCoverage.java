package gridss;

import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialCoverageAnnotator;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;

import java.util.List;
import java.util.concurrent.ExecutorService;

public class AnnotateReferenceCoverage extends VcfTransformCommandLineProgram {
	/**
	 * Defensive programming safety margin around expected window size
	 */
	private final int WINDOW_SIZE_SAFETY_MARGIN = 100000;
	@Override
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		ProcessingContext context = getContext();
		List<SAMEvidenceSource> sources = getSamEvidenceSources();
		int windowSize = SAMEvidenceSource.maximumWindowSize(context, sources, null);
		return new SequentialCoverageAnnotator<VariantContextDirectedEvidence>(context, sources, calls, 2 * windowSize + WINDOW_SIZE_SAFETY_MARGIN, threadpool);
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateReferenceCoverage().instanceMain(argv));
    }
}
