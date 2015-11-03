package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

public class SAMRecordAssemblyEvidenceReadPairIterator extends AbstractIterator<SAMRecordAssemblyEvidence> implements CloseableIterator<SAMRecordAssemblyEvidence>, TrackedBuffer {
	private final ProcessingContext processContext;
	private final AssemblyEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final Iterator<SAMRecord> rit;
	private final SequentialAssemblyReadPairFactory factory;
	private final boolean includeRemote;
	public SAMRecordAssemblyEvidenceReadPairIterator(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			Iterator<SAMRecord> it,
			Iterator<SAMRecord> mateIt,
			boolean includeRemote) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
		this.rit = mateIt;
		this.factory = new SequentialAssemblyReadPairFactory(Iterators.peekingIterator(mateIt));
		this.includeRemote = includeRemote;
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		while (it.hasNext()) {
			SAMRecord record = it.next();
			if (record.getFirstOfPairFlag() || includeRemote) {
				SAMRecordAssemblyEvidence evidence = factory.createAssembly(record, processContext, source);
				if (evidence != null && !evidence.isReferenceAssembly()) {
					return evidence;
				}
			}
		}
		return endOfData();
	}
	@Override
	public void close() {
		CloserUtil.close(it);
		CloserUtil.close(rit);
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