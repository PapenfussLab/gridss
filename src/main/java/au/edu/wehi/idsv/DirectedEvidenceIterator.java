package au.edu.wehi.idsv;

import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Queue;

public class DirectedEvidenceIterator implements CloseableIterator<DirectedEvidence>, PeekingIterator<DirectedEvidence> {
	private static final Log log = Log.getInstance(DirectedEvidenceIterator.class);
	private final SAMEvidenceSource source;
	private Iterator<SAMRecord> it;
	private Queue<DirectedEvidence> buffer = new ArrayDeque<>();
	private int minIndelSize;
	/**
	 * Creates a new iterator
	 * @param it input records
	 */
	public DirectedEvidenceIterator(Iterator<SAMRecord> it, SAMEvidenceSource source, int minIndelSize) {
		if (source == null) {
			throw new IllegalArgumentException("source cannot be null");
		}
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
		if (record == null || record.getReadUnmappedFlag() || record.getMappingQuality() < source.getContext().getConfig().minMapq) {
			return;
		}
		buffer.addAll(SingleReadEvidence.createEvidence(source, minIndelSize, record));
		if (!record.getSupplementaryAlignmentFlag()) {
			if (record.getReadPairedFlag()) {
				ReadPairConcordanceCalculator rpcc = source.getReadPairConcordanceCalculator();
				if (rpcc == null) {
					String msg = String.format("Sanity check failure: null ReadPairConcordanceCalculator for an input file where %s has read paired flag set.", record.getReadName());
					log.error(msg);
					throw new RuntimeException(msg);
				}
				if (!record.getReadUnmappedFlag()
						&& !rpcc.isConcordant(record)) {
					NonReferenceReadPair nrrp = NonReferenceReadPair.create(source, record);
					if (nrrp != null) {
						buffer.add(nrrp);
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
