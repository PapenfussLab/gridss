package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.List;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

public class ReadPairEvidenceIterator extends AbstractIterator<NonReferenceReadPair> implements CloseableIterator<NonReferenceReadPair>, TrackedBuffer {
	private final SAMEvidenceSource source;
	private final SequentialNonReferenceReadPairFactory factory;
	private final Iterator<SAMRecord> iterator;
	private final Iterator<SAMRecord> mateIterator;
	public ReadPairEvidenceIterator(SAMEvidenceSource source, Iterator<SAMRecord> it, Iterator<SAMRecord> mateIt) {
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
			NonReferenceReadPair pair = factory.createNonReferenceReadPair(record, source);
			if (pair != null) {
				return pair;
			}
		}
		return endOfData();
	}
	@Override
	public void setTrackedBufferContext(String context) {
		factory.setTrackedBufferContext(context);
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return factory.currentTrackedBufferSizes();
	}
}
