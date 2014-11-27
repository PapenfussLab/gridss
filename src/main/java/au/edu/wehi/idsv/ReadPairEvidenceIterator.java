package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

public class ReadPairEvidenceIterator extends AbstractIterator<NonReferenceReadPair> implements CloseableIterator<NonReferenceReadPair> {
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final SequentialNonReferenceReadPairFactory factory;
	private final Iterator<SAMRecord> iterator;
	private final Iterator<SAMRecord> mateIterator;
	public ReadPairEvidenceIterator(ProcessingContext processContext, SAMEvidenceSource source, Iterator<SAMRecord> it, Iterator<SAMRecord> mateIt) {
		this.processContext = processContext;
		this.source = source;
		this.iterator = it;
		this.mateIterator = mateIt;
		this.factory = new SequentialNonReferenceReadPairFactory(Iterators.peekingIterator(mateIterator));
	}
	@Override
	public void close() {
		CloserUtil.close(iterator);
		CloserUtil.close(mateIterator);
	}
	@Override
	protected NonReferenceReadPair computeNext() {
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			if (NonReferenceReadPair.meetsLocalEvidenceCritera(processContext.getReadPairParameters(), source, record)) {
				NonReferenceReadPair pair = factory.createNonReferenceReadPair(record, source, processContext);
				if (pair != null && pair.meetsEvidenceCritera(processContext.getReadPairParameters())) {
					return pair;
				}
			}
		}
		return endOfData();
	}
}
