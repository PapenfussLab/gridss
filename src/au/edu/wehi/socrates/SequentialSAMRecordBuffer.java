package au.edu.wehi.socrates;

import com.google.common.collect.PeekingIterator;

import net.sf.samtools.SAMRecord;

public class SequentialSAMRecordBuffer {
	private PeekingIterator<SAMRecord> reads;
	private PeekingIterator<SAMRecord> mates;
	private int referenceIndex;
	public SequentialSAMRecordBuffer(int referenceIndex, PeekingIterator<SAMRecord> reads, PeekingIterator<SAMRecord> mates) {
		this.referenceIndex = referenceIndex;
		this.reads = reads;
		this.mates = mates;
	}
	public SAMRecord getNextEndingSoftClippedRead(int position) {
		// TODO Auto-generated method stub
		return null;
	}

	public NonReferenceReadPair getNextNonReferenceForwardReadPair(int position) {
		// TODO Auto-generated method stub
		return null;
	}

	public SAMRecord getNextStartingSoftClippedRead(int i) {
		// TODO Auto-generated method stub
		return null;
	}

	public NonReferenceReadPair getNextNonReferenceBackwardReadPair(int i) {
		// TODO Auto-generated method stub
		return null;
	}
	public boolean hasReads() {
		// TODO Auto-generated method stub
		return false;
	}
}
