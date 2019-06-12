package au.edu.wehi.validation;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.IntervalUtil;
import htsjdk.samtools.*;
import htsjdk.samtools.SamReaderFactory.Option;

import java.io.File;

public class LongReadSupportFinder {
	private final SamReader reader;
	private final SAMSequenceDictionary dict;
	public LongReadSupportFinder(File indexedBam) {
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(Option.CACHE_FILE_BASED_INDEXES)
				.validationStringency(ValidationStringency.LENIENT);
		this.reader = factory.open(indexedBam);
		this.dict = reader.getFileHeader().getSequenceDictionary();
	}
	/**
	 * evaluate the long read support for the given deletion event
	 * @param chr
	 * @param start1
	 * @param end1
	 * @param start2
	 * @param end2
	 * @param length length of deletion
	 * @param softClipMargin error margin allowed for soft clip location
	 * @param spanningWindowMargin
	 * @param minDeletionSize
	 * @return
	 */
	public LongReadSupportLevel evaluateDeletion(String chr, int start1, int end1, int start2, int end2,
			int softClipMargin,
			int minSoftClipLength,
			int spanningWindowMargin,
			int minDeletionSize) {
		chr = translateReference(chr);
		if (chr == null) return null;
		SAMRecordIterator it = null;
		try {
			it = reader.query(chr, start1, end1, false);
			LongReadSupportLevel support = new LongReadSupportLevel();
			while (it.hasNext()) {
				SAMRecord r = it.next();
				if (SAMRecordUtil.getEndSoftClipLength(r) >= minSoftClipLength && IntervalUtil.overlapsClosed(r.getAlignmentEnd(), r.getAlignmentEnd(), start1 - softClipMargin, end1 + softClipMargin)) {
					support.startClipLocations.add(r.getAlignmentEnd());
				}
				int deletions = countDeletions(r, minDeletionSize, start1, end2);
				if (deletions > 0) {
					support.spanningAlignments.add(deletions);
				}
			}
			it.close();
			it = null;
			it = reader.query(chr, start2, end2, false);
			while (it.hasNext()) {
				SAMRecord r = it.next();
				if (SAMRecordUtil.getStartSoftClipLength(r) >= minSoftClipLength && IntervalUtil.overlapsClosed(r.getAlignmentStart(), r.getAlignmentStart(), start2 - softClipMargin, end2 + softClipMargin)) {
					support.endClipLocations.add(r.getAlignmentStart());
				}
			}
			it.close();
			it = null;
			return support;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		} finally {
			if (it != null) {
				it.close();
			}
		}
	}
	private String translateReference(String chr) {
		if (dict.getSequenceIndex(chr) >= 0) return chr;
		if (dict.getSequenceIndex("chr" + chr) >= 0) return "chr" + chr;
		if (dict.getSequenceIndex(chr.replace("chr",  "")) >= 0) return chr.replace("chr",  "");
		if (dict.getSequenceIndex("Chr" + chr) >= 0) return "Chr" + chr;
		if (dict.getSequenceIndex(chr.replace("Chr",  "")) >= 0) return chr.replace("Chr",  "");
		return null;
	}
	private int countDeletions(SAMRecord r, int minDeletionSize, int startPosition, int endPosition) {
		int position = r.getAlignmentStart();
		int count = 0;
		for (CigarElement ce : r.getCigar().getCigarElements()) {
			if (position >= endPosition) return count;
			if (position >= startPosition && position <= endPosition) {
				if (ce.getOperator() == CigarOperator.D && ce.getLength() > minDeletionSize && position + ce.getLength() < endPosition) {
					// deletion contained within the interval in question
					count += ce.getLength();
				}
			}
			if (ce.getOperator().consumesReferenceBases()) {
				position += ce.getLength();
			}
		}
		return count;
	}
}
