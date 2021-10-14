package au.edu.wehi.idsv;

/**
 * Direction of breakpoint on the given chromosome.
 * <p>The positional relationship between the direct and breakpoint position
 * matches that used by VCF.</p>
 * @author Daniel Cameron
 *
 */
public enum BreakendDirection {
	/**
	 * The breakpoint includes reference bases at and before the breakpoint position
	 * but not after.
	 * 
	 * AAAA.
	 *    ^
	 */
	Forward('f', '+'),
	/**
	 * The breakpoint includes reference bases after but not at the breakpoint position.
	 * 
	 * .AAAA
	 * ^
	 */
	Backward('b', '-');
	private BreakendDirection(char internalChar, char bedChar) {
		this.internalChar = internalChar;
		this.bedChar = bedChar;
	}
	private final char internalChar;
	private final char bedChar;
	public char toChar() { return internalChar; }
	public char toBedChar() { return bedChar; }
	/**
	 * @return Opposite direction of the given breakend direction 
	 */
	public BreakendDirection reverse() { return this == Forward ? Backward : Forward; }
	public static BreakendDirection fromChar(char c) {
		for (BreakendDirection dir : values()) {
			if (dir.toChar() == c) return dir;
		}
		return null;
	}
}
