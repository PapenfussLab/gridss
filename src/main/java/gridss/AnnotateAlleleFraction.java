package gridss;

import au.edu.wehi.idsv.AlleleFractionAnnotator;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import com.google.common.collect.Iterators;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;

import java.util.Iterator;
import java.util.concurrent.ExecutorService;

public class AnnotateAlleleFraction extends VcfTransformCommandLineProgram {
	private final AlleleFractionAnnotator ann;
	public AnnotateAlleleFraction(AlleleFractionAnnotator ann) {
		this.ann = ann;
	}
	@Override
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		Iterator<VariantContextDirectedEvidence> it = Iterators.transform(calls, e -> ann.annotate(e));
		return new AutoClosingIterator<>(it, calls);
	}
}
