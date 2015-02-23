package au.edu.wehi.idsv.alignment;

public interface Aligner {
	/**
	 * Performs Smith-Waterman alignment of the given sequence against the given reference 
	 * @param seq sequence to align
	 * @param ref reference sequence
	 * @return Alignment of sequence relative to reference
	 */
	public Alignment align_smith_waterman(byte[] seq, byte[] ref);
}
