package au.edu.wehi.idsv;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMRecord;

import java.util.Collection;


public interface DirectedEvidence {
	/**
	 * Phred-scaled quality score of breakend
	 * @return
	 */
	float getBreakendQual();
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
	public byte[] getAnchorSequence();
	public byte[] getAnchorQuality();
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
	 * Unique identifier for the source DNA fragments.
	 * @return distinct read names of supporting reads
	 */
	Collection<String> getOriginatingFragmentID(int category);
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
	 * Indicates whether the breakend sequence is exact
	 * @return true if the breakend is known exactly, false otherwise
	 */
	boolean isBreakendExact();
	/**
	 * Strand bias of evidence.
	 * 1 indicates that breakend bases would be aligned to the positive strand if the reference was changed to the variant allele.
	 * 0 indicates that breakend bases would be aligned to the negative strand if the reference was changed to the variant allele.
	 * For read alignments, this corresponds the strand of the read alignment. A SR that is aligned to the positive strand will have a strand bias of 1.
	 * For read pairs, this is the opposite strand to the local anchoring read alignment.
	 * A DP that is aligned to the positive strand will have a strand bias of 0. 
	 * @return Strand bias value between 0 and 1.  
	 */
	double getStrandBias();
	/**
	 * Number of reads/read pairs this evidence composed of.
	 * @return
	 */
	int constituentReads();
	/**
	 * Name of assembly that this read is associated with.
	 * For assembly evidence, this is the name of the assembly itself.
	 * @return Name of associated assembly.
	 */
	String getAssociatedAssemblyName();
	SAMRecord getUnderlyingSAMRecord();
	Ordering<DirectedEvidence> ByEndStart = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return BreakendSummary.ByEndStart.nullsFirst().compare(arg0.getBreakendSummary(), arg1.getBreakendSummary());
		}
	};
	Ordering<DirectedEvidence> ByStartEnd = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return BreakendSummary.ByStartEnd.nullsFirst().compare(arg0.getBreakendSummary(), arg1.getBreakendSummary());
		}
	};
	Ordering<DirectedEvidence> ByBreakendQual = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
			return Doubles.compare(arg0.getBreakendQual(), arg1.getBreakendQual());
		}
	};
	Ordering<DirectedEvidence> ByBreakendQualDesc = ByBreakendQual.reverse();
	Ordering<DirectedEvidence> ByEvidenceID = Ordering.natural().onResultOf(de -> de.getEvidenceID());
}
