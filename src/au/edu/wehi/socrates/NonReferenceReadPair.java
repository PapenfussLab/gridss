package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

/**
 * A read pair that does not support the reference sequence. This can be an OEA, or DP read pair.
 * @author Daniel Cameron
 *
 */
public class NonReferenceReadPair implements DirectedEvidence {
	private final SAMRecord local;
	private final SAMRecord remote;
	private final BreakpointLocation location;
	public NonReferenceReadPair(SAMRecord local, SAMRecord remote, int maxfragmentSize) {
		assert local != null;
		assert remote != null;
		assert (local.getReadName() == null && remote.getReadName() == null) || local.getReadName().equals(remote.getReadName());
		assert !local.getReadUnmappedFlag();
		assert !local.getProperPairFlag();
		this.local = local;
		this.remote = remote;
		this.location = calculateBreakpointLocation(local, remote, maxfragmentSize);
	}
	private static float calculateQualityOEA(SAMRecord local, SAMRecord remote) {
		return local.getMappingQuality();
	}
	private static float calculateQualityDP(SAMRecord local, SAMRecord remote) {
		return Math.min(local.getMappingQuality(), remote.getMappingQuality());
	}
	/**
	 * Calculates the local breakpoint location
	 * @param local local read
	 * @param remote remote read
	 * @param maxfragmentSize maximum fragment size
	 * @return local {@link BreakpointLocation}, without quality information
	 */
	private static BreakpointLocation calculateLocalBreakpointLocation(SAMRecord local, SAMRecord remote, int maxfragmentSize) {
		BreakpointDirection direction = getBreakpointDirection(local);
		int positionClosestToBreakpoint;
		int intervalDirection;
		// adds back in any soft-clipped bases
		int intervalExtendedReadDueToLocalClipping;
		int intervalReducedDueToRemoteMapping = 1;
		if (direction == BreakpointDirection.Forward) {
			positionClosestToBreakpoint = local.getAlignmentEnd();
			intervalDirection = 1;
			intervalExtendedReadDueToLocalClipping = local.getUnclippedEnd() - local.getAlignmentEnd();
		} else {
			positionClosestToBreakpoint = local.getAlignmentStart();
			intervalDirection = -1;
			intervalExtendedReadDueToLocalClipping = local.getAlignmentStart() - local.getUnclippedStart();
		}
		if (remote != null && !remote.getReadUnmappedFlag()) {
			intervalReducedDueToRemoteMapping = remote.getReadLength();
			// add back in any soft-clipped bases
			if (getBreakpointDirection(remote) == BreakpointDirection.Forward) {
				intervalReducedDueToRemoteMapping -= remote.getUnclippedEnd() - remote.getAlignmentEnd();
			} else {
				intervalReducedDueToRemoteMapping -= remote.getAlignmentStart() - remote.getUnclippedStart();
			}
		}
		int intervalWidth = maxfragmentSize - local.getReadLength() + intervalExtendedReadDueToLocalClipping - intervalReducedDueToRemoteMapping;
		return new BreakpointLocation(local.getReferenceIndex(), direction,
				Math.min(positionClosestToBreakpoint, positionClosestToBreakpoint + intervalWidth * intervalDirection),
				Math.max(positionClosestToBreakpoint, positionClosestToBreakpoint + intervalWidth * intervalDirection),
				calculateQualityOEA(local, remote));
	}
	private static BreakpointLocation calculateBreakpointLocation(SAMRecord local, SAMRecord remote, int maxfragmentSize) {
		if (remote == null || remote.getReadUnmappedFlag()) {
			return calculateLocalBreakpointLocation(local, remote, maxfragmentSize);
		} else {
			return new BreakpointInterval(
					calculateLocalBreakpointLocation(local, remote, maxfragmentSize),
					calculateLocalBreakpointLocation(remote, local, maxfragmentSize),
					calculateQualityDP(local, remote));
		}
	}
	/**
	 * Breakpoint direction the read pair supports relative to the given mapped read.
	 * <p>A forward breakpoint direction indicates that this read pairs supports a breakpoint
	 * after the final mapped based of the locally mapped read.</p>
	 * <p>A backward breakpoint direction indicates that this read pairs supports a breakpoint
	 * before the alignment start position of the locally mapped read.</p>
	 * <p>This method assumes an Illumina FR read pair library preparation.</p>
	 * @return breakpoint direction this read supports
	 */
	public static BreakpointDirection getBreakpointDirection(SAMRecord read) {
		return read.getReadNegativeStrandFlag() ? BreakpointDirection.Backward : BreakpointDirection.Forward;
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
	public int getRemoteReferenceIndex() {
		if (remote == null || remote.getReadUnmappedFlag()) return -1;
		return remote.getReferenceIndex();
	}
	@Override
	public String getEvidenceID() {
		return local.getReadName();
	}
	@Override
	public BreakpointLocation getBreakpointLocation() {
		return location;
	}
}
