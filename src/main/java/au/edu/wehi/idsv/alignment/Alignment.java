package au.edu.wehi.idsv.alignment;

/**
 * Alignment of a sequence relative to a reference
 * @author Daniel Cameron
 *
 */
public class Alignment {
	private final int startPosition;
	private final String cigar;
	/**
	 * Sequence pair alignment outcome
	 * @param startPosition Zero-based start position of sequence relative to reference
	 * @param cigar sequence alignment CIGAR
	 */
	public Alignment(int startPosition, String cigar) {
		this.startPosition = startPosition;
		this.cigar = cigar;
	}
	/**
	 * Zero-based start position of sequence relative to reference
	 * @return
	 */
	public int getStartPosition() {
		return startPosition;
	}
	/**
	 * sequence alignment CIGAR
	 * @return CIGAR string of sequence alignment
	 */
	public String getCigar() {
		return cigar;
	}
}
