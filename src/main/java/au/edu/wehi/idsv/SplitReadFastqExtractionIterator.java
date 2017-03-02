package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Queue;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

/**
 * Extracts split read fastq records from the given alignments 
 *
 */
public class SplitReadFastqExtractionIterator implements Iterator<FastqRecord> {
	private final Queue<FastqRecord> buffer = new ArrayDeque<>();
	private final Iterator<SAMRecord> it;
	private final boolean isSplit;
	private final int minSoftClipLength;
	private final float minClipQuality;
	private final boolean processSecondaryAlignments;
	public SplitReadFastqExtractionIterator(
			Iterator<SAMRecord> it,
			boolean isSplit,
			int minSoftClipLength,
			float minClipQuality,
			boolean processSecondaryAlignments) {
		this.it = it;
		this.isSplit = isSplit;
		this.minSoftClipLength = minSoftClipLength;
		this.minClipQuality = minClipQuality;
		this.processSecondaryAlignments = processSecondaryAlignments;
	}

	private void ensureBuffer() {
		while (buffer.isEmpty() && it.hasNext()) {
			SAMRecord r = it.next();
			if (r.getReadUnmappedFlag()) continue;
			// Logic for extending an existing SA alignment not yet complete. Need to:
			// - only realign bases not in any existing SA alignment
			// - update all SA record (requires queryname sorted input file)
			if (r.getAttribute(SAMTag.SA.name()) != null) continue;
			if (r.getSupplementaryAlignmentFlag()) continue;
			if (r.getNotPrimaryAlignmentFlag() && !processSecondaryAlignments) {
				continue;
			}
			for (FastqRecord fqr : SplitReadIdentificationHelper.getSplitReadRealignments(r, isSplit)) {
				if (fqr.length() < minSoftClipLength) continue;
				if (averageBaseQuality(fqr) < minClipQuality) continue;
				buffer.add(fqr);
			}
		}
	}
	private static double averageBaseQuality(FastqRecord fqr) {
		long sum = 0;
		for (byte v : SAMUtils.fastqToPhred(fqr.getBaseQualityString())) {
			sum += v;
		}
		return (double)sum / fqr.getBaseQualityString().length();
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
