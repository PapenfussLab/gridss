package au.edu.wehi.idsv.picard;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

import java.io.IOException;

public class SynchronousReferenceLookupAdapter implements ReferenceLookup {
	private final ReferenceSequenceFile underlying;
	public SynchronousReferenceLookupAdapter(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
	}
	
	@Override
	public synchronized SAMSequenceDictionary getSequenceDictionary() {
		return underlying.getSequenceDictionary();
	}
	
	@Override
	public synchronized ReferenceSequence nextSequence() {
		return underlying.nextSequence();
	}

	@Override
	public synchronized void reset() {
		underlying.reset();
	}

	@Override
	public synchronized boolean isIndexed() {
		return underlying.isIndexed();
	}

	@Override
	public synchronized ReferenceSequence getSequence(String contig) {
		return underlying.getSequence(contig);
	}

	@Override
	public synchronized ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
		return underlying.getSubsequenceAt(contig, start, stop);
	}

	@Override
	public synchronized void close() throws IOException {
		underlying.close();
	}

	@Override
	public synchronized byte getBase(int referenceIndex, int position) {
		return getSubsequenceAt(getSequenceDictionary().getSequence(referenceIndex).getSequenceName(), position, position).getBases()[0];
	}
}
