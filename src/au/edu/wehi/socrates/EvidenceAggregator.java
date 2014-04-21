package au.edu.wehi.socrates;

public class EvidenceClusterProcessor {
	private int referenceIndexFirst;
	private int referenceIndexSecond;
	/**
	 * Calculates reciprocal evidence for breakpoints 
	 * @param referenceIndexLow first chromosome
	 * @param referenceIndexHigh second chromosome. Must have equal or higher BAM reference index
	 */
	public EvidenceClusterProcessor(int referenceIndexFirst, int referenceIndexSecond) {
		this.referenceIndexFirst = referenceIndexFirst;
		this.referenceIndexSecond = referenceIndexSecond;
	}
	// Import all evidence into giant lookup
	
	// lookup contains:
	// supporting evidence: (start, end, qual)
	// breakpoint evidence: (start, end, qual, targetstart, targetend)
	public void load(Iterator<SAMRecord> in) {
	}
}
