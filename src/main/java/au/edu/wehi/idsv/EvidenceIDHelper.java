package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

public class EvidenceIDHelper {
	// can't use /1 and /2 since bwa strips it from the read name during alignment
	// so realignment will fail
	private static final char SEPERATOR = '#'; 
	/**
	 * Gets a string that is unique for the given alignment of the read.
	 * Note: the identifier is only unique with respect to the alignment of the segment/read
	 * not the entire template/read pair.
	 * @param record record
	 * @return alignment-unique identifier for the given SAMRecord
	 */
	public static String getAlignmentUniqueName(SAMRecord record) {
		return buildAlignmentUniqueName(record).toString();
	}
	/**
	 * Extracts an identifier unique for the given read alignment
	 * @param evidenceId evidenceID of evidence from read with the given alignment
	 * @return read alignment unique identifier
	 */
	public static String extractAlignmentUniqueName(String evidenceId) {
		// strip off final seperator
		return evidenceId.substring(0, evidenceId.lastIndexOf(SEPERATOR));
	}
	/**
	 * Extracts an identifier unique for the given read
	 * @param evidenceId evidenceID of evidence from read
	 * @return read alignment unique identifier
	 */
	public static String extractSegmentUniqueName(String evidenceId) {
		return stripSeperators(evidenceId, 5);
	}
	private static String stripSeperators(String evidenceId, int seperatorsToStrip) {
		// strip off all 5 seperators
		int end = evidenceId.length() - 1;
		int seperatorCount = 0;
		while (end >= 0 && seperatorCount < seperatorsToStrip) {
			if (evidenceId.charAt(end) == SEPERATOR) {
				seperatorCount++;
			}
			end--;
		}
		end++;
		if (seperatorCount != seperatorsToStrip) {
			throw new IllegalArgumentException(evidenceId + " is not a valid evidenceID");
		}
		return evidenceId.substring(0, end);
	}
	/**
	 * Extracts the read this evidence was sourced from
	 * @param evidenceId evidenceID of evidence from read
	 * @return read alignment unique identifier
	 */
	public static String extractReadName(String evidenceId) {
		return stripSeperators(evidenceId, 6);
	}
	private static StringBuilder buildSegmentUniqueName(SAMRecord record) {
		StringBuilder sb = new StringBuilder(record.getReadName());
		sb.append(SEPERATOR);
		sb.append(SAMRecordUtil.getSegmentIndex(record));
		return sb;
	}
	private static StringBuilder buildAlignmentUniqueName(SAMRecord record) {
		StringBuilder sb = buildSegmentUniqueName(record);
		if (!record.getReadUnmappedFlag()) {
			sb.append(SEPERATOR);
			sb.append(record.getReferenceName());
			sb.append(SEPERATOR);
			sb.append(record.getAlignmentStart());
			sb.append(SEPERATOR);
			sb.append(record.getReadNegativeStrandFlag() ? '-' : '+');
			sb.append(SEPERATOR);
			sb.append(record.getCigarString());
		} else {
			sb.append(SEPERATOR);
			sb.append(SEPERATOR);
			sb.append(SEPERATOR);
			sb.append(SEPERATOR);
		}
		return sb;
	}
	public static String getEvidenceID(NonReferenceReadPair e) {
		StringBuilder sb = buildAlignmentUniqueName(e.getLocalledMappedRead());
		sb.append(SEPERATOR);
		sb.append("rp");
		// not technically required if only considering 2 segment templates (ie read pairs)
		// but useful for consistency
		if (e.getBreakendSummary() != null) {
			sb.append(e.getBreakendSummary().direction.toChar());
		}
		return sb.toString();
	}
	public static String getEvidenceID(SoftClipEvidence e) {
		StringBuilder sb = buildAlignmentUniqueName(e.getSAMRecord());
		sb.append(SEPERATOR);
		sb.append("sc");
		sb.append(e.getBreakendSummary().direction.toChar());
		return sb.toString();
	}
	public static String getEvidenceID(SplitReadEvidence e) {
		StringBuilder sb = buildAlignmentUniqueName(e.getSAMRecord());
		sb.append(SEPERATOR);
		sb.append("sr");
		sb.append(e.getBreakendSummary().direction.toChar());
		return sb.toString();
	}
	public static String getEvidenceID(IndelEvidence e) {
		StringBuilder sb = buildAlignmentUniqueName(e.getSAMRecord());
		sb.append(SEPERATOR);
		sb.append(e.getIndelCigarOffset());
		sb.append('i');
		sb.append(e.getBreakendSummary().direction.toChar());
		return sb.toString();
	}
}
