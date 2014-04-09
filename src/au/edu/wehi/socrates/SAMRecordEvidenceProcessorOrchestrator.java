package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

public class SAMRecordEvidenceProcessorOrchestrator {
	private SequentialSAMRecordBuffer itbuffer;
	private int maxFragmentSize;
	public SAMRecordEvidenceProcessorOrchestrator(
		SequentialSAMRecordBuffer itbuffer,
		int maxFragmentSize) {
		this.itbuffer = itbuffer;
		this.maxFragmentSize = maxFragmentSize;
	}
	public void process() {
		int position = 1;
		while (itbuffer.hasReads()) {
			// using VCF breakpoint position definition
			// Forward evidence
			for (SAMRecord r = itbuffer.getNextEndingSoftClippedRead(position); r != null; r = itbuffer.getNextEndingSoftClippedRead(position)) {
				// add to forward working set
			}
			for (NonReferenceReadPair r = itbuffer.getNextNonReferenceForwardReadPair(position); r != null; r = itbuffer.getNextNonReferenceForwardReadPair(position)) {
				// add to forward working set
			}
			// Backward evidence
			for (SAMRecord r = itbuffer.getNextStartingSoftClippedRead(position + 1); r != null; r = itbuffer.getNextStartingSoftClippedRead(position + 1)) {
				// add to backward working set
			}
			for (NonReferenceReadPair r = itbuffer.getNextNonReferenceBackwardReadPair(position + maxFragmentSize); r != null; r = itbuffer.getNextNonReferenceBackwardReadPair(position + maxFragmentSize)) {
				// add to backward working set
			}
			// process genomic position
			// advance genomic position to next callable position
			position++;
			// flush expired reads from working sets
		}
	}
}
