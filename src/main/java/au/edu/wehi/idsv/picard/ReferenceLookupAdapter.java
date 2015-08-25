package au.edu.wehi.idsv.picard;

import java.io.IOException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ReferenceLookupAdapter implements ReferenceLookup {
	private final ReferenceSequenceFile underlying;
	public ReferenceLookupAdapter(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
	}
	
	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return underlying.getSequenceDictionary();
	}
	
	@Override
	public ReferenceSequence nextSequence() {
		return underlying.nextSequence();
	}

	@Override
	public void reset() {
		underlying.reset();
	}

	@Override
	public boolean isIndexed() {
		return underlying.isIndexed();
	}

	@Override
	public ReferenceSequence getSequence(String contig) {
		return underlying.getSequence(contig);
	}

	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
		return underlying.getSubsequenceAt(contig, start, stop);
	}

	@Override
	public void close() throws IOException {
		underlying.close();
	}

	@Override
	public byte getBase(int referenceIndex, int position) {
		return getSubsequenceAt(getSequenceDictionary().getSequence(referenceIndex).getSequenceName(), position, position).getBases()[0];
	}
}
