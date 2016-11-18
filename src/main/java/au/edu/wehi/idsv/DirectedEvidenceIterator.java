package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

public class DirectedEvidenceIterator implements CloseableIterator<DirectedEvidence>, PeekingIterator<DirectedEvidence> {
	private final SAMEvidenceSource source;
	private Iterator<SAMRecord> it;
	private Queue<DirectedEvidence> buffer = new ArrayDeque<>();
	/**
	 * Creates a new iterator
	 * @param it input records
	 */
	public DirectedEvidenceIterator(Iterator<SAMRecord> it, SAMEvidenceSource source) {
		this.source = source;
		this.it = it;
	}
	private void ensureBuffer() {
		if (it == null) return;
		if (!buffer.isEmpty()) return;
		while (it != null && it.hasNext() && buffer.isEmpty()) {
			addToBuffer(it.next());
		}
	}
	private void addToBuffer(SAMRecord record) {
		if (record.getReadUnmappedFlag()) return;
		List<SplitReadEvidence> srlist = SplitReadEvidence.create(source, record);
		buffer.addAll(srlist);
		
		// only add soft clip if there isn't a split read
		boolean hasForwardSR = false;
		boolean hasBackwardSR = false;
		for (SplitReadEvidence sre : srlist) {
			switch (sre.getBreakendSummary().direction) {
			case Forward:
				hasForwardSR = true;
				break;
			case Backward:
				hasBackwardSR = true;
				break;
			}
		}
		if (!hasForwardSR && SAMRecordUtil.getEndSoftClipLength(record) > 0) {
			buffer.add(SoftClipEvidence.create(source, BreakendDirection.Forward, record));
		}
		if (!hasBackwardSR && SAMRecordUtil.getStartSoftClipLength(record) > 0) {
			buffer.add(SoftClipEvidence.create(source, BreakendDirection.Backward, record));
		}
		buffer.addAll(IndelEvidence.create(source, record));
		NonReferenceReadPair nrrp = NonReferenceReadPair.create(source, record);
		if (nrrp != null) {
			buffer.add(nrrp);
		}
	}

	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty();
	}

	@Override
	public DirectedEvidence next() {
		if (!hasNext()) throw new NoSuchElementException();
		return buffer.poll();
	}

	@Override
	public DirectedEvidence peek() {
		if (!hasNext()) throw new NoSuchElementException();
		return buffer.peek();
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	@Override
	public void close() {
		CloserUtil.close(it);
		it = null;
	}

}
