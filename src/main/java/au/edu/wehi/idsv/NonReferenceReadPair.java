package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

/**
 * A read pair that does not support the reference sequence. This can be an OEA, or DP read pair.
 * @author Daniel Cameron
 *
 */
public abstract class NonReferenceReadPair implements DirectedEvidence {
	private final SAMRecord local;
	private final SAMRecord remote;
	private final BreakendSummary location;
	private final SAMEvidenceSource source;
	protected NonReferenceReadPair(SAMRecord local, SAMRecord remote, SAMEvidenceSource source) {
		if (local == null) throw new IllegalArgumentException("local is null");
		if (remote == null) throw new IllegalArgumentException("remote is null");
		if (!StringUtils.equals(local.getReadName(), remote.getReadName())) throw new IllegalArgumentException(String.format("Paired reads %s and %s have differing read names", local.getReadName(), remote.getReadName()));
		if (local.getReadUnmappedFlag()) throw new IllegalArgumentException("local must be mapped");
		//if (local.getProperPairFlag()) throw new IllegalArgumentException(String.format("Read %s is flagged as part of a proper pair", local.getReadName()));
		//if (remote.getProperPairFlag()) throw new IllegalArgumentException(String.format("Read %s is flagged as part of a proper pair", remote.getReadName()));
		if (source.getMaxConcordantFragmentSize() < local.getReadLength()) throw new IllegalArgumentException(String.format("Sanity check failure: read pair %s contains read of length %d when maximum fragment size is %d", local.getReadName(), local.getReadLength(), source.getMaxConcordantFragmentSize()));
		this.local = local;
		this.remote = remote;
		this.location = calculateBreakendSummary(local, remote, source.getMaxConcordantFragmentSize());
		this.source = source;
	}
	public static NonReferenceReadPair create(SAMRecord local, SAMRecord remote, SAMEvidenceSource source) {
		if (remote == null || remote.getReadUnmappedFlag()) {
			return new UnmappedMateReadPair(local, remote, source);
		} else {
			return new DiscordantReadPair(local, remote, source);
		}
	}
	public boolean meetsEvidenceCritera(ReadPairParameters rp) {
		return location != null
				&& meetsLocalEvidenceCritera(rp, source, local)
				&& meetsRemoteEvidenceCritera(rp, source, remote)
				&& !SAMRecordUtil.isDovetailing(local,  remote);
	}
	public static boolean meetsLocalEvidenceCritera(ReadPairParameters rp, SAMEvidenceSource source, SAMRecord local) {
		return !local.getReadUnmappedFlag()
				&& local.getReadPairedFlag()
				&& local.getMappingQuality() >= rp.minLocalMapq
				&& !SAMRecordUtil.estimatedReadsOverlap(local)
				&& !source.getReadPairConcordanceCalculator().isConcordant(local);
	}
	public static boolean meetsRemoteEvidenceCritera(ReadPairParameters rp, SAMEvidenceSource source, SAMRecord remote) {
		return remote.getReadPairedFlag()
				&& !remote.getMateUnmappedFlag()
				&& ((remote.getReadUnmappedFlag() && !remote.getMateUnmappedFlag()) // OEA
					|| !(source.getReadPairConcordanceCalculator().isConcordant(remote) || SAMRecordUtil.estimatedReadsOverlap(remote))); // DP
	}
	/**
	 * Calculates the local breakpoint location
	 * @param local local read
	 * @param remote remote read
	 * @param maxfragmentSize maximum fragment size
	 * @return local {@link BreakendSummary}, without quality information
	 */
	private static BreakendSummary calculateLocalBreakendSummary(SAMRecord local, SAMRecord remote, int maxfragmentSize) {
		BreakendDirection direction = getBreakendDirection(local);
		int positionClosestToBreakpoint;
		int intervalDirection;
		// adds back in any soft-clipped bases
		int intervalExtendedReadDueToLocalClipping;
		int intervalReducedDueToRemoteMapping = 1;
		if (direction == BreakendDirection.Forward) {
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
			if (getBreakendDirection(remote) == BreakendDirection.Forward) {
				intervalReducedDueToRemoteMapping -= remote.getUnclippedEnd() - remote.getAlignmentEnd();
			} else {
				intervalReducedDueToRemoteMapping -= remote.getAlignmentStart() - remote.getUnclippedStart();
			}
		}
		int intervalWidth = maxfragmentSize - local.getReadLength() + intervalExtendedReadDueToLocalClipping - intervalReducedDueToRemoteMapping;
		intervalWidth = Math.min(intervalWidth, pairSeparation(local, remote));
		if (intervalWidth < 0) return null;
		return new BreakendSummary(local.getReferenceIndex(), direction,
				Math.min(positionClosestToBreakpoint, positionClosestToBreakpoint + intervalWidth * intervalDirection),
				Math.max(positionClosestToBreakpoint, positionClosestToBreakpoint + intervalWidth * intervalDirection));
	}
	/**
	 * Determines the separation between discordant reads
	 * Determines the number of unsequenced bases in the fragment
	 * @param local
	 * @param remote
	 * @return number possible breakpoints between the read pair mapped in the expected orientation,
	 *  or Integer.MAX_VALUE if placement is not as expected
	 */
	private static int pairSeparation(SAMRecord local, SAMRecord remote) {
		if (local.getReadUnmappedFlag() || remote.getReadUnmappedFlag()) return Integer.MAX_VALUE;
		if (local.getReferenceIndex() != remote.getReferenceIndex()) return Integer.MAX_VALUE;
		// Assuming FR orientation
		if (local.getReadNegativeStrandFlag() == remote.getReadNegativeStrandFlag()) return Integer.MAX_VALUE;
				// <--local
		if ((local.getReadNegativeStrandFlag() && local.getAlignmentStart() > remote.getAlignmentStart())
				// local--> 
				|| (!local.getReadNegativeStrandFlag() && local.getAlignmentStart() < remote.getAlignmentStart())) {
			// only problem with this pair is the fragment size is unexpected
			return Math.max(local.getAlignmentStart(), remote.getAlignmentStart()) - Math.min(local.getAlignmentEnd(), remote.getAlignmentEnd()) - 1;
		}
		return Integer.MAX_VALUE;
	}
	private static BreakendSummary calculateBreakendSummary(SAMRecord local, SAMRecord remote, int maxfragmentSize) {
		if (remote == null || remote.getReadUnmappedFlag()) {
			return calculateLocalBreakendSummary(local, remote, maxfragmentSize);
		} else {
			// Discordant because the pairs overlap = no SV evidence
			int separation = pairSeparation(local, remote);
			if (separation < 0) return null;
			BreakendSummary bsLocal = calculateLocalBreakendSummary(local, remote, maxfragmentSize);
			BreakendSummary bsRemote = calculateLocalBreakendSummary(remote, local, maxfragmentSize);
			if (bsLocal == null || bsRemote == null) return null;
			return new BreakpointSummary(bsLocal, bsRemote);
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
	public static BreakendDirection getBreakendDirection(SAMRecord read) {
		return read.getReadNegativeStrandFlag() ? BreakendDirection.Backward : BreakendDirection.Forward;
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
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	@Override
	public SAMEvidenceSource getEvidenceSource() {
		return source;
	}
	@Override
	public int getLocalMapq() {
		return local.getMappingQuality();
	}
	@Override
	public int getLocalBaseLength() {
		return local.getReadLength();
	}
	@Override
	public int getLocalBaseCount() {
		return local.getReadLength();
	}
	@Override
	public int getLocalMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(local);
	}
	@Override
	public int getLocalTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(local);
	}
	@Override
	public byte[] getBreakendSequence() {
		return null;
	}
	@Override
	public byte[] getBreakendQuality() {
		return null;
	}
}
