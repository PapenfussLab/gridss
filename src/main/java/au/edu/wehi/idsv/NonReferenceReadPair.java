package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamPairUtil.PairOrientation;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

/**
 * A read pair that does not support the reference sequence. This can be an OEA, or DP read pair.
 * @author Daniel Cameron
 *
 */
public abstract class NonReferenceReadPair implements DirectedEvidence {
	//private static final Log log = Log.getInstance(NonReferenceReadPair.class);
	private final SAMRecord local;
	private final SAMRecord remote;
	private final BreakendSummary location;
	private final SAMEvidenceSource source;
	private String evidenceID = null;
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
		this.location = calculateBreakendSummary(local, meetsAnchorCriteria(source, remote) ? remote : null, source);
		this.source = source;
	}
	/**
	 * Creates a read pair not supporting the reference allele from the given reads
	 * @param local read mapped to local side of putative breakpoint
	 * @param remote read mapped to remote side of putative breakpoint
	 * @param source source
	 * @return Read pair evidence if reads support putative breakpoint, null otherwise
	 */
	public static NonReferenceReadPair create(SAMRecord local, SAMRecord remote, SAMEvidenceSource source) {
		if (local == null || remote == null) return null;
		if (SAMRecordUtil.entropy(local) < source.getContext().getConfig().minAnchorShannonEntropy) return null;
		if (SAMRecordUtil.entropy(remote) < source.getContext().getConfig().minAnchorShannonEntropy) return null;
		if (!meetsAnchorCriteria(source, local)) return null;
		if (SAMRecordUtil.isDovetailing(local,  remote, PairOrientation.FR, source.getContext().getConfig().dovetailMargin)) return null; 
		// should only need to check for adapters in OEA as DP with adapter should be dovetailing
		if (source.getContext().getConfig().adapters.containsAdapter(local)) return null;
		if (source.getContext().getConfig().adapters.containsAdapter(remote)) return null;
		assert(source.getReadPairConcordanceCalculator().isConcordant(local, remote) == source.getReadPairConcordanceCalculator().isConcordant(remote, local));
		if (source.getReadPairConcordanceCalculator().isConcordant(local, remote)) return null;
		if (pairSeparation(local, remote, PairOrientation.FR) < 0) return null; // discordant because the pairs overlap = no SV evidence
		if (SAMRecordUtil.areSameRead(local, remote)) {
			//log.debug(String.format("Filtering self-matching read %s at position %s:%d", local.getReadName(), source.getContext().getDictionary().getSequence(local.getReferenceIndex()).getSequenceName(), local.getAlignmentStart()));
			return null;
		}
		NonReferenceReadPair rp = null;
		if (!meetsAnchorCriteria(source, remote)) {
			rp = new UnmappedMateReadPair(local, remote, source);
		} else {
			rp = new DiscordantReadPair(local, remote, source);
		}
		if (rp.location == null) {
			// all seemed good but we found no SV support
			// Eg: read length >= expected max fragment size 
			rp = null;
		}
		return rp;
	}
	public static boolean meetsAnchorCriteria(SAMEvidenceSource source, SAMRecord read) {
		return read.getReadPairedFlag()
				&& !read.getReadUnmappedFlag()
				&& !read.getReadFailsVendorQualityCheckFlag()
				&& read.getMappingQuality() >= source.getContext().getConfig().minMapq
				&& read.getMappingQuality() <= source.getContext().getConfig().maxMapq
				&& !SAMRecordUtil.estimatedReadsOverlap(read, PairOrientation.FR, source.getMetrics().getIdsvMetrics().MAX_READ_LENGTH - source.getMetrics().getMaxSoftClipLength())
				&& !source.getReadPairConcordanceCalculator().isConcordant(read);
	}
	public static boolean meetsRemoteCriteria(SAMEvidenceSource source, SAMRecord read) {
		return read.getReadPairedFlag()
				&& !read.getMateUnmappedFlag() // other side has to be an anchor
				&& !read.getReadFailsVendorQualityCheckFlag()
				&& (meetsDPRemoteCriteria(source, read) || meetsOEARemoteCriteria(source, read))
				&& !source.getContext().getConfig().adapters.containsAdapter(read);
	}
	private static boolean meetsDPRemoteCriteria(SAMEvidenceSource source, SAMRecord read) {
		return read.getReadPairedFlag()
			&& !read.getReadUnmappedFlag()
			&& !read.getReadFailsVendorQualityCheckFlag()
			&& !SAMRecordUtil.estimatedReadsOverlap(read, PairOrientation.FR, source.getMetrics().getIdsvMetrics().MAX_READ_LENGTH - source.getMetrics().getMaxSoftClipLength())
			&& !source.getReadPairConcordanceCalculator().isConcordant(read);
	}
	private static boolean meetsOEARemoteCriteria(SAMEvidenceSource source, SAMRecord read) {
		return read.getReadUnmappedFlag();
	}
	/**
	 * Calculates the local breakpoint location
	 * @param local local read
	 * @param remote remote read
	 * @param maxfragmentSize maximum fragment size
	 * @return local {@link BreakendSummary}, without quality information
	 */
	private static BreakendSummary calculateLocalBreakendSummary(SAMRecord local, SAMRecord remote, int maxfragmentSize, SAMSequenceDictionary dictionary) {
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
		intervalWidth = Math.min(intervalWidth, pairSeparation(local, remote, PairOrientation.FR));
		if (intervalWidth < 0) return null;
		return new BreakendSummary(local.getReferenceIndex(), direction,
				Math.max(1, Math.min(positionClosestToBreakpoint, positionClosestToBreakpoint + intervalWidth * intervalDirection)),
				Math.min(dictionary.getSequence(local.getReferenceIndex()).getSequenceLength(), Math.max(positionClosestToBreakpoint, positionClosestToBreakpoint + intervalWidth * intervalDirection)));
	}
	/**
	 * Determines the separation between discordant reads
	 * Determines the number of unsequenced bases in the fragment
	 * @param local
	 * @param remote
	 * @param expectedOrientation 
	 * @return number possible breakpoints between the read pair mapped in the expected orientation,
	 *  or Integer.MAX_VALUE if placement is not as expected
	 */
	private static int pairSeparation(SAMRecord local, SAMRecord remote, PairOrientation expectedOrientation) {
		if (remote == null || local.getReadUnmappedFlag() || remote.getReadUnmappedFlag()) return Integer.MAX_VALUE;
		if (local.getReferenceIndex() != remote.getReferenceIndex()) return Integer.MAX_VALUE;
		// Assuming FR orientation
		if (local.getReadNegativeStrandFlag() == remote.getReadNegativeStrandFlag()) return Integer.MAX_VALUE;
				// <--local
		if ((local.getReadNegativeStrandFlag() && local.getAlignmentStart() >= remote.getAlignmentStart())
				// local--> 
				|| (!local.getReadNegativeStrandFlag() && local.getAlignmentStart() <= remote.getAlignmentStart())) {
			// only problem with this pair is the fragment size is unexpected
			return Math.max(local.getAlignmentStart(), remote.getAlignmentStart()) - Math.min(local.getAlignmentEnd(), remote.getAlignmentEnd()) - 1;
		}
		return Integer.MAX_VALUE;
	}
	private static BreakendSummary calculateBreakendSummary(SAMRecord local, SAMRecord remote, SAMEvidenceSource source) {
		int maxFragmentSize = source.getMaxConcordantFragmentSize();
		SAMSequenceDictionary dictionary = source.getContext().getDictionary();
		BreakendSummary bsLocal = calculateLocalBreakendSummary(local, remote, maxFragmentSize, dictionary);
		if (remote == null || remote.getReadUnmappedFlag()) {
			return bsLocal;
		} else {
			BreakendSummary bsRemote = calculateLocalBreakendSummary(remote, local, maxFragmentSize, dictionary);
			if (bsRemote == null && bsLocal == null) return null;
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
	public static String getEvidenceID(SAMRecord record) {
		return record.getReadName() + (record.getFirstOfPairFlag() ? SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX : SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX);
	}	
	@Override
	public String getEvidenceID() {
		if (evidenceID == null) {
			evidenceID = getEvidenceID(local);
		}
		return evidenceID;
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
	public byte[] getBreakendSequence() {
		return null;
	}
	@Override
	public byte[] getBreakendQuality() {
		return null;
	}
	@Override
	public byte[] getAnchorSequence() {
		return null;
	}
	@Override
	public byte[] getAnchorQuality() {
		return null;
	}
	@Override
	public boolean isBreakendExact() {
		return false;
	}
	/**
	 * Determines whether the reads are on strands compatible with their expected orientation.
	 * Alignment details are ignored
	 * @return true if the strands are as expected, false if the strands are unexpected. 
	 */
	public boolean onExpectedStrand() {
		PairOrientation po = PairOrientation.FR;
		if (source.getMetrics() != null && source.getMetrics().getInsertSizeMetrics() != null && source.getMetrics().getInsertSizeMetrics().PAIR_ORIENTATION != null) {
			po = source.getMetrics().getInsertSizeMetrics().PAIR_ORIENTATION;
		}
		switch (po) {
			case FR:
			case RF:
				return getLocalledMappedRead().getReadNegativeStrandFlag() != getNonReferenceRead().getReadNegativeStrandFlag();
			case TANDEM:
				return getLocalledMappedRead().getReadNegativeStrandFlag() == getNonReferenceRead().getReadNegativeStrandFlag();
		}
		throw new IllegalStateException("Unknown read pair orientation");
	}
}