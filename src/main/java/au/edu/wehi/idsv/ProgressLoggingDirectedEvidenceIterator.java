package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;

import java.util.Iterator;

public class ProgressLoggingDirectedEvidenceIterator<T extends DirectedEvidence> implements CloseableIterator<T> {
	private final ProgressLoggerInterface logger;
	private final Iterator<T> iterator;
	private final GenomicProcessingContext processContext;
	public ProgressLoggingDirectedEvidenceIterator(GenomicProcessingContext processContext, Iterator<T> iterator, ProgressLoggerInterface logger) {
		this.iterator = iterator;
		this.logger = logger;
		this.processContext = processContext;
	}
	@Override
	public boolean hasNext() {
		return iterator.hasNext();
	}
	@Override
	public T next() {
		T n = iterator.next();
		BreakendSummary bs = n.getBreakendSummary();
		if (bs != null) {
			logger.record(processContext.getDictionary().getSequence(bs.referenceIndex).getSequenceName(), bs.start);
		}
		return n;
	}
	@Override
	public void close() {
		CloserUtil.close(iterator);
	}
}
