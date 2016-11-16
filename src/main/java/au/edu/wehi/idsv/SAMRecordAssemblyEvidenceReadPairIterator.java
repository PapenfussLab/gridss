package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

public class SAMRecordAssemblyEvidenceReadPairIterator implements CloseableIterator<SAMRecordAssemblyEvidence>, TrackedBuffer {
	private final ProcessingContext processContext;
	private final AssemblyEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final Iterator<SAMRecord> rit;
	private final SequentialAssemblyReadPairFactory factory;
	private final boolean includeRemote;
	private final boolean includeSpannedIndels;
	private Queue<SAMRecordAssemblyEvidence> buffer = new ArrayDeque<SAMRecordAssemblyEvidence>(); 
	public SAMRecordAssemblyEvidenceReadPairIterator(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			Iterator<SAMRecord> it,
			Iterator<SAMRecord> mateIt,
			boolean includeRemote,
			boolean includeSpannedIndels) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
		this.rit = mateIt;
		this.factory = new SequentialAssemblyReadPairFactory(Iterators.peekingIterator(mateIt));
		this.includeRemote = includeRemote;
		this.includeSpannedIndels = includeSpannedIndels;
	}
	@Override
	public void close() {
		CloserUtil.close(it);
		CloserUtil.close(rit);
		buffer.clear();
	}
	@Override
	public void setTrackedBufferContext(String context) {
		factory.setTrackedBufferContext(context);
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return factory.currentTrackedBufferSizes();
	}
	private void ensureBuffer() {
		while (buffer.isEmpty() && it.hasNext()) {
			SAMRecord record = it.next();
			if (record.getFirstOfPairFlag() || includeRemote) {
				SAMRecordAssemblyEvidence evidence = factory.createAssembly(record, processContext, source);
				if (evidence != null) {
					if (includeSpannedIndels && !(evidence instanceof RemoteEvidence)) {
						buffer.addAll(evidence.getSpannedIndels());
					}
					if (evidence.getBreakendSummary() != null) {
						buffer.add(evidence);
					}
				}
			}
		}
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty(); 
	}
	@Override
	public SAMRecordAssemblyEvidence next() {
		ensureBuffer();
		return buffer.poll();
	}
}