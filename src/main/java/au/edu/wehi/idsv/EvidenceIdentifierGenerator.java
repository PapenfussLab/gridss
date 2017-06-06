package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public interface EvidenceIdentifierGenerator { 
	/**
	 * Gets a string that is unique for the given alignment of the given read.
	 * @param record record
	 * @return alignment-unique identifier for the given SAMRecord
	 */
	String getAlignmentUniqueName(SAMRecord record);
	/**
	 * Gets a string that is unique for the read segment. That is, 
	 * for paired-end sequencing, this ensures that the first a second
	 * reads of a read pair have unique identifiers
	 * @param record record
	 * @return segment-unique identifier for the given SAMRecord
	 */
	String getSegmentUniqueName(SAMRecord record);
	/**
	 * Extracts an identifier unique for the given read alignment
	 * @param evidenceId evidenceID of evidence from read with the given alignment
	 * @return read alignment unique identifier
	 */
	String extractAlignmentUniqueName(String evidenceId);
	/**
	 * Extracts an identifier unique for the given read
	 * @param evidenceId evidenceID of evidence from read
	 * @return read alignment unique identifier
	 */
	String extractSegmentUniqueName(String evidenceId);
	/**
	 * Extracts the read this evidence was sourced from
	 * @param evidenceId evidenceID of evidence from read
	 * @return read alignment unique identifier
	 */
	//public static String extractReadName(String evidenceId) {
	//	return stripSeperators(evidenceId, 6);
	//}
	String getEvidenceID(NonReferenceReadPair e);
	String getEvidenceID(SoftClipEvidence e);
	String getEvidenceID(SplitReadEvidence e);
	String getEvidenceID(IndelEvidence e);
}
