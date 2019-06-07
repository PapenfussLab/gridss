package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Queue;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.PeekingIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

public class DirectedEvidenceIterator implements CloseableIterator<DirectedEvidence>, PeekingIterator<DirectedEvidence> {
	private final SAMEvidenceSource source;
	private Iterator<SAMRecord> it;
	private Queue<DirectedEvidence> buffer = new ArrayDeque<>();
	private int minIndelSize;
	/**
	 * Creates a new iterator
	 * @param it input records
	 */
	public DirectedEvidenceIterator(Iterator<SAMRecord> it, SAMEvidenceSource source, int minIndelSize) {
		this.source = source;
		this.it = it;
		this.minIndelSize = minIndelSize;
	}
	private void ensureBuffer() {
		if (it == null) return;
		if (!buffer.isEmpty()) return;
		while (it != null && it.hasNext() && buffer.isEmpty()) {
			addToBuffer(it.next());
		}
	}
	private void addToBuffer(SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getMappingQuality() < source.getContext().getConfig().minMapq) {
			return;
		}
		buffer.addAll(SingleReadEvidence.createEvidence(source, minIndelSize, record));
		if (!record.getSupplementaryAlignmentFlag()) {
			if (NonReferenceReadPair.meetsAnchorCriteria(source, record)) {
				NonReferenceReadPair nrrp = NonReferenceReadPair.create(source, record);
				if (nrrp != null) {
					buffer.add(nrrp);
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
