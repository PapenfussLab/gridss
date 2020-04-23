package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Extracts split read fastq records from the given alignments 
 *
 */
public class SplitReadFastqExtractor {
	private final boolean isSplit;
	private final int minSoftClipLength;
	private final float minClipQuality;
	private final boolean processSecondaryAlignments;
	private final boolean realignExistingSplitReads;
	private final boolean realignEntireRecord;
	private final boolean realignAnchoringBases;
	private final EvidenceIdentifierGenerator eidgen;
	public SplitReadFastqExtractor(
			boolean isSplit,
			int minSoftClipLength,
			float minClipQuality,
			boolean processSecondaryAlignments,
			boolean realignExistingSplitReads,
			boolean realignEntireRecord,
			boolean realignAnchoringBases,
			EvidenceIdentifierGenerator eidgen) {
		this.isSplit = isSplit;
		this.minSoftClipLength = minSoftClipLength;
		this.minClipQuality = minClipQuality;
		this.processSecondaryAlignments = processSecondaryAlignments;
		this.realignExistingSplitReads = realignExistingSplitReads;
		this.realignEntireRecord = realignEntireRecord;
		this.realignAnchoringBases = realignAnchoringBases;
		this.eidgen = eidgen;
		if (realignEntireRecord && realignAnchoringBases) {
			throw new IllegalArgumentException("Cannot realign anchoring bases if realigning entire read");
		}
	}
	public List<FastqRecord> extract(SAMRecord r) {
		List<FastqRecord> list = new ArrayList<>(2);
		if (r.getReadUnmappedFlag()) return list;
		if (!SAMRecordUtil.isSoftClipLengthAtLeast(r, minSoftClipLength)) return list;
		if (realignEntireRecord && !isSplit && !AssemblyAttributes.isUnanchored(r)) {
			list.add(SplitReadHelper.getFullRealignment(r, eidgen));
			return list;
		}
		if (isSplit && SplitReadHelper.isAnchoringBasesRecord(r)) {
			// don't chain realigning the anchoring bases alignments record - it's a once off
			return list;
		}
		// Logic for extending an existing SA alignment not yet complete. Need to:
		// - only realign bases not in any existing SA alignment
		// - update all SA record (requires queryname sorted input file)
		// Note that realignments may be split reads, but we currently only consider the primary alignment
		if (!realignExistingSplitReads && !isSplit && r.getAttribute(SAMTag.SA.name()) != null) return list;
		if (r.getSupplementaryAlignmentFlag()) return list;
		if (r.isSecondaryAlignment() && !processSecondaryAlignments) return list;
		for (FastqRecord fqr : SplitReadHelper.getSplitReadRealignments(r, isSplit, eidgen)) {
			if (fqr.getReadLength() < minSoftClipLength) continue;
			if (averageBaseQuality(fqr) < minClipQuality) continue;
			list.add(fqr);
		}
		if (realignAnchoringBases && !isSplit && !AssemblyAttributes.isUnanchored(r) && list.size() > 0) {
			// only emit the anchoring bases:
			// - on the first pass where the anchoring bases are actually the read anchoring bases
			// - when we actually have anchoring bases
			// - when we have something that we might actually realign (TODO: is this the best approach?)
			list.add(SplitReadHelper.getAnchoringBases(r, eidgen));
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
