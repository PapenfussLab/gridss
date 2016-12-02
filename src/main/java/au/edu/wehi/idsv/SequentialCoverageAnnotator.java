package au.edu.wehi.idsv;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;

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
 * @author Daniel Cameron
 *
 */
public class SequentialCoverageAnnotator<T extends VariantContextDirectedEvidence> implements CloseableIterator<T> {
	private static final Log log = Log.getInstance(SequentialCoverageAnnotator.class);
	private final ProcessingContext context;
	private final List<ReferenceCoverageLookup> reference;
	private final Iterator<T> it;
	private final List<Closeable> toclose = new ArrayList<>();
	public SequentialCoverageAnnotator(ProcessingContext context, List<SAMEvidenceSource> sources, Iterator<T> it) {
		int windowSize = SAMEvidenceSource.maximumWindowSize(context, sources, null);
		List<List<SAMEvidenceSource>> byCategory = Lists.newArrayList();
		while (byCategory.size() < context.getCategoryCount()) {
			byCategory.add(new ArrayList<SAMEvidenceSource>());
		}
		for (SAMEvidenceSource s : sources) {
			byCategory.get(s.getSourceCategory()).add(s);
		}
		List<ReferenceCoverageLookup> coverage = Lists.newArrayList();
		for (int i = 0; i < byCategory.size(); i++) {
			List<ReferenceCoverageLookup> lookup = Lists.newArrayList();
			for (SAMEvidenceSource s : byCategory.get(i)) {
				// one read-ahead thread per input file
				SamReader reader = SamReaderFactory.makeDefault().open(s.getFile());
				SAMRecordIterator rawIterator = reader.iterator();
				rawIterator.assertSorted(SortOrder.coordinate);
				CloseableIterator<SAMRecord> sit = new AsyncBufferedIterator<SAMRecord>(rawIterator, s.getFile().getName() + "-Coverage");
				toclose.add(reader);
				toclose.add(rawIterator);
				toclose.add(sit);
				sit = new ProgressLoggingSAMRecordIterator(sit, new ProgressLogger(log));
				SequentialReferenceCoverageLookup sourceLookup = new SequentialReferenceCoverageLookup(sit, s.getMetrics().getIdsvMetrics(), s.getReadPairConcordanceCalculator(), windowSize);
				context.registerBuffer(s.getFile().getName(), sourceLookup);
				lookup.add(sourceLookup);
			}
			coverage.add(new AggregateReferenceCoverageLookup(lookup));
		}
		this.context = context;
		this.reference = coverage;
		this.it = it;
		
	}
	public SequentialCoverageAnnotator(
			ProcessingContext context,
			Iterator<T> it,
			List<ReferenceCoverageLookup> reference) {
		this.it = it;
		this.context = context;
		this.reference = reference;
	}
	@SuppressWarnings("unchecked")
	public T annotate(T variant) {
		BreakendSummary loc = variant.getBreakendSummary();
		int referenceIndex = loc.referenceIndex;
		int offset = loc.direction == BreakendDirection.Forward ? 0 : -1;
		int[] reads = new int[reference.size()];
		int[] spans = new int[reference.size()];
		for (int i = 0; i < reference.size(); i++) {
			ReferenceCoverageLookup r = reference.get(i);
			if (r != null) {
				reads[i] = IntStream.range(loc.start + offset, loc.end + 1 + offset)
						.map(p -> r.readsSupportingNoBreakendAfter(referenceIndex, p))
						.min().getAsInt();
				spans[i] = IntStream.range(loc.start + offset, loc.end + 1 + offset)
						.map(p -> r.readPairsSupportingNoBreakendAfter(referenceIndex, p))
						.min().getAsInt();
			}
		}
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(context, variant);
		builder.referenceReads(reads);
		builder.referenceSpanningPairs(spans);
		return (T)builder.make();
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
