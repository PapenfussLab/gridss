package au.edu.wehi.idsv;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Soft clip split read in which the realigned soft clip maps to a local coordinate
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedRemoteSoftClipEvidence extends RealignedSoftClipEvidence implements RemoteEvidence {
	private final RealignedSoftClipEvidence local;
	public RealignedRemoteSoftClipEvidence(RealignedSoftClipEvidence evidence) {
		super(
				evidence.getEvidenceSource(),
				getRemoteDirection(evidence.getBreakendSummary().direction, evidence.getSAMRecord(), evidence.getRealignedSAMRecord()),
				createLocal(evidence, evidence.getBreakendSummary().direction, evidence.getSAMRecord(), evidence.getRealignedSAMRecord()),
				createRemote(evidence, evidence.getBreakendSummary().direction, evidence.getSAMRecord(), evidence.getRealignedSAMRecord()));
		this.local = evidence;
		assert(getBreakendSummary().equals(local.getBreakendSummary().remoteBreakpoint()));
	}
	public RealignedSoftClipEvidence asLocal() {
		return local;
	}
	@Override
	public String getEvidenceID() {
		return getSAMRecord().getReadName();
	}
	private static String getRemoteEvidenceID(RealignedSoftClipEvidence evidence) {
		return "R" + evidence.getEvidenceID();
	}
	/**
	 * Creates a new SAMRecord equivalent to a realigned record had to aligner mapped to the realigned location
	 * @param direction
	 * @param local
	 * @param realigned
	 * @return
	 */
	private static SAMRecord createRemote(RealignedSoftClipEvidence evidence, BreakendDirection direction, SAMRecord local, SAMRecord realigned) {
		// strip out the soft clip to convert into a realigned record
		SAMRecord newRealigned = SAMRecordUtil.clone(local);
		newRealigned.setReadName(getRemoteEvidenceID(evidence));
		List<CigarElement> cigar = Lists.newArrayList(newRealigned.getCigar().getCigarElements());
		int readLen = local.getReadLength();
		byte[] bases = local.getReadBases();
		byte[] quals = local.getBaseQualities();
		int startOffset, endOffset;
		// untemplated sequence soft clip gets moved over to this new realigned side
		int untemplated = getUntemplatedSequenceLength(direction, local, realigned);
		int newLength = local.getReadLength() - realigned.getReadLength() + untemplated;
		if (direction == BreakendDirection.Forward) {
			cigar.remove(cigar.size() - 1);
			if (untemplated != 0) {
				cigar.add(new CigarElement(untemplated, CigarOperator.SOFT_CLIP));
			}
			startOffset = 0;
			endOffset = newLength;
		} else {
			cigar.remove(0);
			if (untemplated != 0) {
				cigar.add(0, new CigarElement(untemplated, CigarOperator.SOFT_CLIP));
			}
			startOffset = readLen - newLength;
			endOffset = readLen;
		}
		newRealigned.setCigar(new Cigar(cigar));
		if (bases != SAMRecord.NULL_SEQUENCE) {
			newRealigned.setReadBases(Arrays.copyOfRange(bases, startOffset, endOffset));
		}
		if (quals != SAMRecord.NULL_QUALS) {
			newRealigned.setBaseQualities(Arrays.copyOfRange(quals, startOffset, endOffset));
		}
		// set breakend orientation
		BreakendDirection remoteDirection = getRemoteDirection(direction, local, realigned);
		newRealigned.setReadNegativeStrandFlag(direction == remoteDirection);
		assert(newRealigned.isValid() == null);
		return newRealigned;
	}
	/**
	 * Creates a new SAMRecord equivalent to a soft clipped record had to aligner mapped to the realigned location
	 * @param evidence 
	 * @param direction
	 * @param local
	 * @param realigned
	 * @return
	 */
	public static SAMRecord createLocal(RealignedSoftClipEvidence evidence, BreakendDirection direction, SAMRecord local, SAMRecord realigned) {
		// fill in the anchor bases as a large soft clip on the realigned record so it acts like a local soft clip record
		SAMRecord newLocal = SAMRecordUtil.clone(realigned);
		newLocal.setReadName(getRemoteEvidenceID(evidence));
		List<CigarElement> cigar = Lists.newArrayList(newLocal.getCigar().getCigarElements());
		newLocal.setReadNegativeStrandFlag(local.getReadNegativeStrandFlag()); // preserve source strand so we can get back to the FASTQ sequence
		byte[] bases = local.getReadBases();
		byte[] quals = local.getBaseQualities();
		if (bases != SAMRecord.NULL_SEQUENCE) {
			if (realigned.getReadNegativeStrandFlag()) {
				bases = bases.clone();
				SequenceUtil.reverseComplement(bases);
			}
			newLocal.setReadBases(bases);
		}
		if (quals != SAMRecord.NULL_QUALS) {
			if (realigned.getReadNegativeStrandFlag()) {
				quals = bases.clone();
				SequenceUtil.reverseComplement(quals);
			}
			newLocal.setBaseQualities(quals);
		}
		int untemplated = getUntemplatedSequenceLength(direction, local, realigned);
		int newScLength = local.getReadLength() - realigned.getReadLength() + untemplated;
		if (getRemoteDirection(direction, local, realigned) == BreakendDirection.Forward) {
			if (untemplated > 0) {
				assert(cigar.get(cigar.size() - 1).getOperator() == CigarOperator.SOFT_CLIP);
				assert(cigar.get(cigar.size() - 1).getLength() == untemplated);
				cigar.remove(cigar.size() - 1);
			}
			cigar.add(new CigarElement(newScLength, CigarOperator.SOFT_CLIP));
		} else {
			if (untemplated > 0) {
				assert(cigar.get(0).getOperator() == CigarOperator.SOFT_CLIP);
				assert(cigar.get(0).getLength() == untemplated);
				cigar.remove(0);
			}
			cigar.add(0, new CigarElement(newScLength, CigarOperator.SOFT_CLIP));
		}
		newLocal.setCigar(new Cigar(cigar));
		assert(newLocal.isValid() == null);
		return newLocal;
	}
	private static BreakendDirection getRemoteDirection(BreakendDirection direction, SAMRecord local, SAMRecord realigned) {
		if (direction == BreakendDirection.Forward) {
			if (realigned.getReadNegativeStrandFlag()) {
				return BreakendDirection.Forward;
			} else {
				return BreakendDirection.Backward;
			}
		} else {
			if (realigned.getReadNegativeStrandFlag()) {
				return BreakendDirection.Backward;
			} else {
				return BreakendDirection.Forward;
			}
		}
	}
	private static int getUntemplatedSequenceLength(BreakendDirection direction, SAMRecord local, SAMRecord realigned) {
		List<CigarElement> ce = realigned.getCigar().getCigarElements();
		if (getRemoteDirection(direction, local, realigned) == BreakendDirection.Forward) {
			if (ce.get(ce.size() - 1).getOperator() == CigarOperator.SOFT_CLIP) {
				return ce.get(ce.size() - 1).getLength();
			}
		} else {
			if (ce.get(0).getOperator() == CigarOperator.SOFT_CLIP) {
				return ce.get(0).getLength();
			}
		}
		return 0;
	}
	@Override
	public String toString() {
		return "R" + super.toString();
	}
	@Override
	public float getBreakpointQual() {
		return local.getBreakpointQual();
	}
	@Override
	public float getBreakendQual() {
		return local.getBreakendQual();
	}
}
