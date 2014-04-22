package au.edu.wehi.socrates;

import org.junit.experimental.max.MaxCore;

import net.sf.samtools.SAMRecord;

/**
 * A read pair that does not support the reference sequence. This can be an OEA, or DP read pair.
 * @author Daniel Cameron
 *
 */
public class NonReferenceReadPair implements DirectedEvidence {
	private final SAMRecord local;
	private final SAMRecord remote;
	private final int fragmentSize;
	public NonReferenceReadPair(SAMRecord local, SAMRecord remote, int fragmentSize) {
		assert local != null;
		assert remote != null;
		assert (local.getReadName() == null && remote.getReadName() == null) || local.getReadName().equals(remote.getReadName());
		assert !local.getReadUnmappedFlag();
		assert !local.getProperPairFlag();
		this.local = local;
		this.remote = remote;
		this.fragmentSize = fragmentSize;
	}
	/**
	 * Mapped read under consideration
	 * @return
	 */
	public SAMRecord getLocalledMappedRead() { return local; }
	/**
	 * Read not supporting the reference placement of the originating fragment 
	 * @return
	 */
	public SAMRecord getNonReferenceRead() { return remote; }
	/**
	 * Breakpoint direction the read pair supports relative to the locally mapped read.
	 * <p>A forward breakpoint direction indicates that this read pairs supports a breakpoint
	 * after the final mapped based of the locally mapped read.</p>
	 * <p>A backward breakpoint direction indicates that this read pairs supports a breakpoint
	 * before the alignment start position of the locally mapped read.</p>
	 * <p>This method assumes an Illumina FR read pair library preparation.</p>
	 * @return breakpoint direction this read pair supports
	 */
	public BreakpointDirection getBreakpointDirection() {
		return local.getReadNegativeStrandFlag() ? BreakpointDirection.Backward : BreakpointDirection.Forward;
	}
	@Override
	public String getEvidenceID() {
		return local.getReadName();
	}
	@Override
	public int getReferenceIndex() {
		return local.getReferenceIndex();
	}
	@Override
	public long getWindowStart() {
		return getBreakpointDirection() == BreakpointDirection.Forward ? local.getAlignmentEnd() : local.getAlignmentStart() - fragmentSize;
	}
	@Override
	public long getWindowEnd() {
		return getBreakpointDirection() == BreakpointDirection.Forward ? local.getAlignmentEnd() + fragmentSize : local.getAlignmentStart();
	}
}
