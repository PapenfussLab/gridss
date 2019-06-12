package au.edu.wehi.idsv;

import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;
import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;
import java.util.Base64;

/**
 * Generates evidence ID using hashes to reduce evidenceID size
 * 
 * The encoding consists three url-safe Base64 encoded byte blocks (6 bits of information per byte):
 * The first block is the segment unique hash (typically 20 byte = 120 bits)
 * The second block is the alignment unique hash for that segment (typically 6 bytes = 36 bits)
 * The final block is the overall evidenceid hash for that alignment (typically 6 bytes = 36 bits)
 * 
 * 
 * @author Daniel Cameron
 *
 */
public class HashedEvidenceIdentifierGenerator implements EvidenceIdentifierGenerator {
	private final StringEvidenceIdentifierGenerator gen = new StringEvidenceIdentifierGenerator();
	private final Base64.Encoder encoder = Base64.getUrlEncoder().withoutPadding();
	//private final Base64.Decoder decoder = Base64.getUrlDecoder();
	private HashFunction hf = Hashing.murmur3_128();
	private final int segmentUniqueBytes;
	private final int alignmentUniqueBytes;
	private final int evidenceidUniqueBytes;
	public HashedEvidenceIdentifierGenerator(int segmentUniqueBytes, int alignmentUniqueBytes, int evidenceidUniqueBytes) {
		this.segmentUniqueBytes = segmentUniqueBytes;
		this.alignmentUniqueBytes = alignmentUniqueBytes;
		this.evidenceidUniqueBytes = evidenceidUniqueBytes;
	}
	public HashedEvidenceIdentifierGenerator() {
		this(20, 6, 6);
	}
	/**
	 * Hashes the given string, returning a string that does not contain any SAM or VCF special characters.
	 * @param s
	 * @return string encoding hash
	 */
	// SAM read name regex: \*|[!-()+-<>-~][!-~]*
	// !"#$%&'()+,-./0123456789:;<>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
	// * = disallowed
	private String hash(String s, int bytes) {
		HashCode hc = hf.hashString(s, StandardCharsets.US_ASCII);
		String encoded = encoder.encodeToString(hc.asBytes());
		String truncated = encoded.substring(0, bytes);
		return truncated;
	}
	@Override
	public String extractAlignmentUniqueName(String evidenceId) {
		return evidenceId.substring(0, segmentUniqueBytes + alignmentUniqueBytes);
	}
	@Override
	public String extractSegmentUniqueName(String evidenceId) {
		return evidenceId.substring(0, segmentUniqueBytes);
	}
	@Override
	public String getAlignmentUniqueName(SAMRecord record) {
		String aun = gen.getAlignmentUniqueName(record);
		return getSegmentUniqueName(record) + hash(aun, alignmentUniqueBytes);
	}
	@Override
	public String getSegmentUniqueName(SAMRecord record) {
		String sun = gen.getSegmentUniqueName(record);
		return hash(sun, segmentUniqueBytes);
	}
	@Override
	public String getEvidenceID(NonReferenceReadPair e) {
		String id = gen.getEvidenceID(e);
		return getAlignmentUniqueName(e.getLocalledMappedRead()) + hash(id, evidenceidUniqueBytes);
	}
	@Override
	public String getEvidenceID(SoftClipEvidence e) {
		String id = gen.getEvidenceID(e);
		return getAlignmentUniqueName(e.getSAMRecord()) + hash(id, evidenceidUniqueBytes);
	}
	@Override
	public String getEvidenceID(SplitReadEvidence e) {
		String id = gen.getEvidenceID(e);
		return getAlignmentUniqueName(e.getSAMRecord()) + hash(id, evidenceidUniqueBytes);
	}
	@Override
	public String getEvidenceID(IndelEvidence e) {
		String id = gen.getEvidenceID(e);
		return getAlignmentUniqueName(e.getSAMRecord()) + hash(id, evidenceidUniqueBytes);
	}
}
