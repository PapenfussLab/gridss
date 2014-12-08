package au.edu.wehi.idsv;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;


public interface DirectedEvidence {
	/**
	 * Location of breakpoints consistent with the given evidence.
	 * If the destination of the breakpoint is known, a @see BreakpointSummary
	 * should be returned. 
	 * @return breakpoint locations implied by this evidence
	 */
	BreakendSummary getBreakendSummary();
	/**
	 * Gets the breakpoint sequence excluding anchored based according to the breakend anchor positive strand.   
	 * @return breakpoint sequence bases, null if breakend sequence is not known
	 */
	public byte[] getBreakendSequence();
	/**
	 * Gets the breakpoint sequence quality
	 * @return 0-based phred-like quality scores, null if breakend sequence is not known
	 */
	public byte[] getBreakendQuality();
	/**
	 * Unique breakpoint identifier.
	 * This identifier is used as the FASTQ sequence identifier
	 * when performing realignment and must be a valid FASTQ
	 * sequence identifier as well as being unique across all
	 * evidence processors.</p>
	 * @return Unique breakpoint identifier string
	 */
	String getEvidenceID();
	/**
	 * Source of this evidence
	 * @return Source providing this evidence
	 */
	EvidenceSource getEvidenceSource();
	/**
	 * Maximum MAPQ of SV-supporting evidence mapped to the reference 
	 * @return
	 */
	int getLocalMapq();
	/**
	 * Length of reference mapped sequence
	 * @return
	 */
	int getLocalBaseLength();
	/**
	 * Maximum base quality of reference mapped bases
	 * @return
	 */
	int getLocalMaxBaseQual();
	/**
	 * Total base quality of reference mapped bases
	 * @return
	 */
	int getLocalTotalBaseQual();
	static final Ordering<DirectedEvidence> ByEndStart = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return BreakendSummary.ByEndStart.compare(arg0.getBreakendSummary(), arg1.getBreakendSummary());
		}
	};
	static final Ordering<DirectedEvidence> ByStartEnd = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return BreakendSummary.ByStartEnd.compare(arg0.getBreakendSummary(), arg1.getBreakendSummary());
		}
	};
	static final Ordering<DirectedEvidence> ByLlr = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return Doubles.compare(Models.llr(arg0), Models.llr(arg1));
		}
	};
	static final Ordering<DirectedEvidence> ByLlrDesc = ByLlr.reverse();
}
