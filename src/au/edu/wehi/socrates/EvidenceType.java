package au.edu.wehi.socrates;

/**
 * Types of structural variation breakpoint evidence
 * @author Daniel Cameron
 *
 */
public enum EvidenceType {
	/**
	 * Read contains soft-clipped bases
	 */
	SoftClip,
	/**
	 * Open-ended anchor read pairs. Only one read in each pair is mapped.
	 */
	OpenEndedAnchor,
	/**
	 * Discordant read pairs
	 */
	DiscordantPair
}
