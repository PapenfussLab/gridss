package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Queue;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.fastq.FastqRecord;

/**
 * Extracts split read fastq records from the given alignments 
 *
 */
public class SplitReadFastqExtractionIterator implements Iterator<FastqRecord> {
	private final Queue<FastqRecord> buffer = new ArrayDeque<>();
	private final Iterator<SAMRecord> it;
	private final boolean isSplit;
	private final boolean processSecondaryAlignments;
	private final int minSoftClipLength;
	
	public SplitReadFastqExtractionIterator(
			Iterator<SAMRecord> it,
			boolean isSplit,
			int minSoftClipLength,
			boolean processSecondaryAlignments) {
		this.it = it;
		this.isSplit = isSplit;
		this.minSoftClipLength = minSoftClipLength;
		this.processSecondaryAlignments = processSecondaryAlignments;
	}

	private void ensureBuffer() {
		while (buffer.isEmpty() && it.hasNext()) {
			SAMRecord r = it.next();
			if (r.getReadUnmappedFlag()) continue;
			// Logic for extending an existing SA alignment not yet complete. Need to:
			// - only realign bases on in any existing SA alignment
			// - update all SA record (requires queryname sorted input file)
			if (r.getAttribute(SAMTag.SA.name()) != null) continue;
			if (r.getSupplementaryAlignmentFlag()) continue;
			if (!processSecondaryAlignments && r.getNotPrimaryAlignmentFlag()) continue;
			for (FastqRecord fqr : SplitReadIdentificationHelper.getSplitReadRealignments(r, isSplit)) {
				if (fqr.length() >= minSoftClipLength) {
					buffer.add(fqr);
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
	public FastqRecord next() {
		ensureBuffer();
		if (buffer.isEmpty()) {
			throw new NoSuchElementException();
		}
		return buffer.poll();
	}
}
