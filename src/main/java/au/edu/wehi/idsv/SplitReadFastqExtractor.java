package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;

/**
 * Extracts split read fastq records from the given alignments 
 *
 */
public class SplitReadFastqExtractor {
	private final boolean isSplit;
	private final int minSoftClipLength;
	private final float minClipQuality;
	private final boolean processSecondaryAlignments;
	private final EvidenceIdentifierGenerator eidgen;
	public SplitReadFastqExtractor(
			boolean isSplit,
			int minSoftClipLength,
			float minClipQuality,
			boolean processSecondaryAlignments,
			EvidenceIdentifierGenerator eidgen) {
		this.isSplit = isSplit;
		this.minSoftClipLength = minSoftClipLength;
		this.minClipQuality = minClipQuality;
		this.processSecondaryAlignments = processSecondaryAlignments;
		this.eidgen = eidgen;
	}
	public List<FastqRecord> extract(SAMRecord r) {
		List<FastqRecord> list = new ArrayList<>(2);
		if (r.getReadUnmappedFlag()) return list;
		// Logic for extending an existing SA alignment not yet complete. Need to:
		// - only realign bases not in any existing SA alignment
		// - update all SA record (requires queryname sorted input file)
		// Note that realignments may be split reads, but we currently only consider the primary alignment
		if (!isSplit && r.getAttribute(SAMTag.SA.name()) != null) return list;
		if (r.getSupplementaryAlignmentFlag()) return list;
		if (r.getNotPrimaryAlignmentFlag() && !processSecondaryAlignments) {
			return list;
		}
		for (FastqRecord fqr : SplitReadIdentificationHelper.getSplitReadRealignments(r, isSplit, eidgen)) {
			if (fqr.length() < minSoftClipLength) continue;
			if (averageBaseQuality(fqr) < minClipQuality) continue;
			list.add(fqr);
		}
		return list;
	}
	private static double averageBaseQuality(FastqRecord fqr) {
		long sum = 0;
		for (byte v : SAMUtils.fastqToPhred(fqr.getBaseQualityString())) {
			sum += v;
		}
		return (double)sum / fqr.getBaseQualityString().length();
	}
}
