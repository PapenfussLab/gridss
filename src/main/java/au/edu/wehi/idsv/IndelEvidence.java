package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import au.edu.wehi.idsv.sam.CigarUtil;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

/**
 * Read alignment based support for a small structural variant
 * @author Daniel Cameron
 *
 */
public class IndelEvidence extends SingleReadEvidence implements DirectedBreakpoint {
	private final List<CigarElement> indel;
	private int indelCigarElementOffset;
	private IndelEvidence remote;
	private IndelEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd,
			int offsetRemoteStart, int offsetRemoteEnd,
			List<CigarElement> indel,
			int indelCigarElementOffset) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, offsetRemoteStart, offsetRemoteEnd, 0, 0);
		this.indel = indel;
		this.indelCigarElementOffset = indelCigarElementOffset;
	}
	public static IndelEvidence create(SAMEvidenceSource source, SAMRecord record, int indelCigarElementOffset) {
		List<CigarElement> cl = record.getCigar().getCigarElements();
		int indelEndOffset = indelCigarElementOffset; // exclusive end offset
		while (cl.get(indelEndOffset).getOperator().isIndelOrSkippedRegion() && indelEndOffset < cl.size()) {
			indelEndOffset++;
		}
		assert(indelEndOffset > indelCigarElementOffset);
		assert(indelCigarElementOffset > 0);
		assert(indelEndOffset < cl.size());
		assert(!cl.get(indelCigarElementOffset - 1).getOperator().isIndelOrSkippedRegion());
		List<CigarElement> pre = cl.subList(0, indelCigarElementOffset);
		List<CigarElement> indel = cl.subList(indelCigarElementOffset, indelEndOffset);
		List<CigarElement> post = cl.subList(indelEndOffset, cl.size());
		int preReadLength = CigarUtil.readLength(pre);
		int indelReadLength = CigarUtil.readLength(indel);
		int postReadLength = CigarUtil.readLength(post);
		int preRefLength = CigarUtil.referenceLength(pre);
		//int indelRefLength = CigarUtil.referenceLength(indel);
		int postRefLength = CigarUtil.referenceLength(post);
		// TODO: adjust bounds based on sequence homology
		BreakpointSummary location = new BreakpointSummary(
				record.getReferenceIndex(), BreakendDirection.Forward, record.getAlignmentStart() + preRefLength - 1,
				record.getReferenceIndex(), BreakendDirection.Backward, record.getAlignmentEnd() - postRefLength + 1);
		int preStartOffset = 0;
		int preEndOffset = preReadLength;
		int unmappedStartOffset = preEndOffset;
		int unmappedEndOffset = unmappedStartOffset + indelReadLength;
		int postStartOffset = unmappedEndOffset;
		int postEndOffset = unmappedEndOffset + postReadLength;
		
		if (!INCLUDE_CLIPPED_ANCHORING_BASES) {
			preStartOffset += CigarUtil.getStartSoftClipLength(pre);
			postEndOffset -= CigarUtil.getEndSoftClipLength(post);
		}
		IndelEvidence left = new IndelEvidence(source, record, location,
				preStartOffset, preEndOffset, unmappedStartOffset, unmappedEndOffset, postStartOffset, postEndOffset, indel, indelCigarElementOffset);
		IndelEvidence right = new IndelEvidence(source, record, location.remoteBreakpoint(),
				postStartOffset, postEndOffset, unmappedStartOffset, unmappedEndOffset, preStartOffset, preEndOffset, indel, indelCigarElementOffset);
		left.remote = right;
		right.remote = left;
		return left;
	}
	public static List<IndelEvidence> create(SAMEvidenceSource source, SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) return Collections.emptyList();
		if (CigarUtil.widthOfImprecision(record.getCigar()) > 0) {
			// not a real indel: this is a placeholder CIGAR for an unanchored breakend assembly
			return Collections.emptyList();
		}
		List<IndelEvidence> list = new ArrayList<IndelEvidence>(4);
		List<CigarElement> cl = record.getCigar().getCigarElements();
		// find all indels (excluding those at start/end)
		int indelStartOffset = 1;
		while (indelStartOffset < cl.size()) {
			if (cl.get(indelStartOffset).getOperator().isIndelOrSkippedRegion() &&
					!cl.get(indelStartOffset - 1).getOperator().isIndelOrSkippedRegion()) {
				int indelEndOffset = indelStartOffset;
				while (indelEndOffset < cl.size() && cl.get(indelEndOffset).getOperator().isIndelOrSkippedRegion()) {
					indelEndOffset++;
				}
				if (indelEndOffset < cl.size()) { 
					IndelEvidence left = create(source, record, indelStartOffset);
					IndelEvidence right = left.asRemote();
					list.add(left);
					list.add(right);
				}
				indelStartOffset = indelEndOffset;
			} else {
				indelStartOffset++;
			}
		}
		return list;
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return (BreakpointSummary)super.getBreakendSummary();
	}
	
	@Override
	public float getBreakendQual() {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return scoreAssembly() / 2;
		}
		CigarElement e = indel.get(0);
		for (int i = 1; i < indel.size(); i++) {
			if (e.getLength() < indel.get(i).getLength()) {
				e = indel.get(i);
			}
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreIndel(getEvidenceSource().getMetrics(), e.getOperator(), e.getLength(), getLocalMapq()) / 2;
	}
	@Override
	public float getBreakpointQual() {
		return getBreakendQual();
	}
	private float scoreAssembly() {
		if (getBreakendSequence().length == 0) return 0;
		AssemblyAttributes attr = new AssemblyAttributes(getSAMRecord());
		int rp = attr.getAssemblySupportCountReadPair();
		double rpq = attr.getAssemblySupportReadPairQualityScore();
		int sc = attr.getAssemblySupportCountSoftClip();
		double scq =  attr.getAssemblySupportSoftClipQualityScore();
		if (source.getContext().getAssemblyParameters().excludeNonSupportingEvidence) {
			rp -= attr.getAssemblyNonSupportingReadPairCount();
			rpq -= attr.getAssemblyNonSupportingReadPairQualityScore();
			sc -= attr.getAssemblyNonSupportingSoftClipCount();
			scq -= attr.getAssemblyNonSupportingSoftClipQualityScore();
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreAssembly(
				rp, rpq,
				sc, scq,
				getLocalMapq(),
				getRemoteMapq());
	}

	@Override
	public int getRemoteMapq() {
		return getLocalMapq();
	}

	@Override
	public IndelEvidence asRemote() {
		assert(remote != null);
		return remote;
	}
	@Override
	public String getRemoteEvidenceID() {
		return remote.getEvidenceID();
	}
	
	@Override
	protected String getUncachedEvidenceID() {
		return EvidenceIDHelper.getEvidenceID(this);
	}
	/**
	 * Identifies which indel in the read this evidence corresponds to.
	 * @return zero-based offset in the read CIGAR operator list of this indel
	 */
	public int getIndelCigarOffset() {
		return indelCigarElementOffset;
	}

}
