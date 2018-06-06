package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Queue;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

/**
 * Extracts split read fastq records from the given alignments 
 *
 */
public class SplitReadFastqExtractionIterator implements Iterator<FastqRecord> {
	private final Queue<FastqRecord> buffer = new ArrayDeque<>();
	private final Iterator<SAMRecord> it;
	private final SplitReadFastqExtractor extractor;
	public SplitReadFastqExtractionIterator(
			Iterator<SAMRecord> it,
			boolean isSplit,
			int minSoftClipLength,
			float minClipQuality,
			boolean processSecondaryAlignments,
			boolean realignExistingSplitReads,
			boolean realignEntireRecord,
			EvidenceIdentifierGenerator eidgen) {
		this(it, new SplitReadFastqExtractor(
				isSplit,
				minSoftClipLength,
				minClipQuality,
				processSecondaryAlignments,
				realignExistingSplitReads,
				realignEntireRecord,
				eidgen));
	}
	
	public SplitReadFastqExtractionIterator(
			Iterator<SAMRecord> it,
			SplitReadFastqExtractor extractor) {
		this.it = it;
		this.extractor = extractor;
	}

	private void ensureBuffer() {
		while (buffer.isEmpty() && it.hasNext()) {
			buffer.addAll(extractor.extract(it.next()));
		}
	}

	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty();
	}

	@Override
	public FastqRecord next() {
		ensureBuffer();
		if (buffer.isEmpty()) {
			throw new NoSuchElementException();
		}
		return buffer.poll();
	}
}
