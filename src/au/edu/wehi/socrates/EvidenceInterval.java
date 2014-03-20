package au.edu.wehi.socrates;
/**
 * Relative interval over which the given evidence provides support for a genomic position 
 * @author DANIEL
 *
 */
public abstract class EvidenceInterval {
	public static EvidenceInterval newEvidenceInterval(EvidenceType type, BreakpointDirection direction, int startOffset, int endOffset) {
		switch (type) {
		}
	}
	private final int startOffset;
	private final int endOffset;
	protected EvidenceInterval(int startOffset, int endOffset) {
		this.startOffset = startOffset;
		this.endOffset = endOffset;
	}
	public int getStartOffset() { return this.startOffset; }
	public int getEndOffset() { return this.endOffset; }
}
