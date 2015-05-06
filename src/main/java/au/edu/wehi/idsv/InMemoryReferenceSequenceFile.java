package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

import java.io.IOException;
import java.util.Arrays;

public class InMemoryReferenceSequenceFile implements ReferenceSequenceFile {
	private SAMSequenceDictionary dictionary;
	private byte[][] sequences;
	private int referenceIndex = 0;
	public InMemoryReferenceSequenceFile(String[] contigNames, byte[][] sequences) {
		dictionary = new SAMSequenceDictionary();
		for (int i = 0; i < contigNames.length; i++) {
			dictionary.addSequence(new SAMSequenceRecord(contigNames[i], sequences[i].length));	
		}
		this.sequences = sequences;
	}
	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return dictionary;
	}
	@Override
	public ReferenceSequence nextSequence() {
		ReferenceSequence rs = null;
		if (referenceIndex <= sequences.length) {
			rs = new ReferenceSequence(dictionary.getSequence(referenceIndex).getSequenceName(), referenceIndex, sequences[referenceIndex]);
			referenceIndex++;
		}
		return rs;
	}
	@Override
	public void reset() {
		referenceIndex = 0;
	}
	@Override
	public boolean isIndexed() {
		return true;
	}
	@Override
	public ReferenceSequence getSequence(String contig) {
		int index = dictionary.getSequence(contig).getSequenceIndex();
		return new ReferenceSequence(contig, index, sequences[index]);
	}
	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
		int index = dictionary.getSequence(contig).getSequenceIndex();
		return new ReferenceSequence(contig, index, Arrays.copyOfRange(sequences[index], (int)start, (int)stop));
	}
	@Override
	public void close() throws IOException {
	}
}
