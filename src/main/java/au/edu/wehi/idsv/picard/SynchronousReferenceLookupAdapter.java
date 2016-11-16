package au.edu.wehi.idsv.picard;

import java.io.IOException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class SynchronousReferenceLookupAdapter implements ReferenceLookup {
	private final ReferenceSequenceFile underlying;
	private boolean isClosed = false;
	public SynchronousReferenceLookupAdapter(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
	}
	
	@Override
	public synchronized SAMSequenceDictionary getSequenceDictionary() {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		return underlying.getSequenceDictionary();
	}
	
	@Override
	public synchronized ReferenceSequence nextSequence() {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		return underlying.nextSequence();
	}

	@Override
	public synchronized void reset() {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		underlying.reset();
	}

	@Override
	public synchronized boolean isIndexed() {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		return underlying.isIndexed();
	}

	@Override
	public synchronized ReferenceSequence getSequence(String contig) {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		return underlying.getSequence(contig);
	}

	@Override
	public synchronized ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		return underlying.getSubsequenceAt(contig, start, stop);
	}

	@Override
	public synchronized void close() throws IOException {
		if (!isClosed) {
			underlying.close();
		}
		isClosed = true;
	}

	@Override
	public synchronized byte getBase(int referenceIndex, int position) {
		if (isClosed) throw new IllegalStateException("Underlying reference has closed");
		return getSubsequenceAt(getSequenceDictionary().getSequence(referenceIndex).getSequenceName(), position, position).getBases()[0];
	}
}
