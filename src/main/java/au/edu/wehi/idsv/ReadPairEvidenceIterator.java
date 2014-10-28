package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

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
			if (SAMRecordUtil.isPartOfNonReferenceReadPair(record)) {
				NonReferenceReadPair pair = factory.createNonReferenceReadPair(record, source);
				if (pair != null && pair.isValid()
						&& processContext.getReadPairParameters().meetsEvidenceCritera(pair)
						// overlapping fragments in the expected orientation do not indicate a SV
						// TODO: what about dovedtailing fragments? Likely due to oversequencing of small fragment but could be SV
						&& !(pair.fragmentSequencesOverlap() && SAMRecordUtil.couldBeProperPair(record, source.getMetrics().getPairOrientation()))) {
					return pair;
				}
			}
		}
		return endOfData();
	}
}
