package au.edu.wehi.idsv;

public interface ExactBreakend extends DirectedEvidence {
	/**
	 * Gets the breakpoint sequence excluding anchor.
	 * @return breakpoint sequence
	 */
	public byte[] getBreakendSequence();

	/**
	 * Gets the breakpoint sequence quality
	 * @return 0-based phred-like quality scores
	 */
	public byte[] getBreakendQuality();

}