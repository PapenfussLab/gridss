package au.edu.wehi.idsv;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;


/**
 * Annotates breakends with reference allele coverage information
 * 
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialCoverageAnnotator<T extends VariantContextDirectedEvidence> implements CloseableIterator<T> {
	private static final Log log = Log.getInstance(SequentialCoverageAnnotator.class);
	private final ProcessingContext context;
	private final List<ReferenceCoverageLookup> reference;
	private final Iterator<T> it;
	private final List<Closeable> toclose = new ArrayList<>();
	private final ExecutorService threadpool;
	public SequentialCoverageAnnotator(ProcessingContext context, List<SAMEvidenceSource> sources, Iterator<T> it, ExecutorService threadpool) {
		this.context = context;
		this.reference = createLookup(context, sources);
		this.it = it;
		this.threadpool = threadpool;
	}
	private List<ReferenceCoverageLookup> createLookup(ProcessingContext context, List<SAMEvidenceSource> sources) {
		int windowSize = SAMEvidenceSource.maximumWindowSize(context, sources, null);
		List<ReferenceCoverageLookup> result = new ArrayList<>();
		for (SAMEvidenceSource ses : sources) {
			assert(ses.getSourceCategory() >= 0);
			assert(ses.getSourceCategory() < context.getCategoryCount());
			// one read-ahead thread per input file
			SamReader reader = SamReaderFactory.makeDefault().open(ses.getFile());
			SAMRecordIterator rawIterator = reader.iterator();
			rawIterator.assertSorted(SortOrder.coordinate);
			CloseableIterator<SAMRecord> sit = new AsyncBufferedIterator<SAMRecord>(rawIterator, ses.getFile().getName() + "-Coverage");
			toclose.add(sit); // close the async iterator first to prevent aysnc reading from a closed stream 
			toclose.add(rawIterator);
			toclose.add(reader);
			sit = new ProgressLoggingSAMRecordIterator(sit, new ProgressLogger(log));
			SequentialReferenceCoverageLookup sourceLookup = new SequentialReferenceCoverageLookup(sit, ses.getMetrics().getIdsvMetrics(), ses.getReadPairConcordanceCalculator(), windowSize, ses.getSourceCategory());
			context.registerBuffer(ses.getFile().getName(), sourceLookup);
			result.add(sourceLookup);
		}
		return result;
	}
	public SequentialCoverageAnnotator(
			ProcessingContext context,
			Iterator<T> it,
			List<ReferenceCoverageLookup> reference,
			ExecutorService threadpool) {
		this.it = it;
		this.context = context;
		this.reference = reference;
		this.threadpool = threadpool;
	}
	private static class CoverageResult {
		public CoverageResult(int reads, int spans) {
			this.readsSupportingNoBreakendAfter = reads;
			this.readPairsSupportingNoBreakendAfter = spans;
		}
		public final int readsSupportingNoBreakendAfter;
		public final int readPairsSupportingNoBreakendAfter;
	}
	private static CoverageResult calculateCoverage(ReferenceCoverageLookup lookup, int referenceIndex, int start, int end) {
		int reads = IntStream.range(start, end)
				.map(p -> lookup.readsSupportingNoBreakendAfter(referenceIndex, p))
				.min().getAsInt();
		int spans = IntStream.range(start, end)
				.map(p -> lookup.readPairsSupportingNoBreakendAfter(referenceIndex, p))
				.min().getAsInt();
		return new CoverageResult(reads, spans);
	}
	@SuppressWarnings("unchecked")
	public T annotate(T variant) {
		BreakendSummary loc = variant.getBreakendSummary();
		int referenceIndex = loc.referenceIndex;
		int offset = loc.direction == BreakendDirection.Forward ? 0 : -1;
		int start = loc.start + offset;
		int end = loc.end + 1 + offset;
		List<Future<CoverageResult>> tasks = new ArrayList<>();
		for (ReferenceCoverageLookup rcl : reference) {
			tasks.add(threadpool.submit(() -> calculateCoverage(rcl, referenceIndex, start, end)));
		}
		try {
			int[] reads = new int[context.getCategoryCount()];
			int[] spans = new int[context.getCategoryCount()];
			for (int i = 0; i < reference.size(); i++) {
				ReferenceCoverageLookup rcl = reference.get(i);
				assert(rcl.getCategory() < context.getCategoryCount());
				reads[rcl.getCategory()] += tasks.get(i).get().readsSupportingNoBreakendAfter;
				spans[rcl.getCategory()] += tasks.get(i).get().readPairsSupportingNoBreakendAfter;
			}
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context, variant);
			builder.referenceReads(reads);
			builder.referenceSpanningPairs(spans);
			return (T)builder.make();
		} catch (ExecutionException | InterruptedException e) {
			throw new RuntimeException(e);
		}
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public T next() {
		return annotate(it.next());
	}
	@Override
	public void close() {
		for (Closeable c : toclose) {
			CloserUtil.close(c);
		}
	}
}
