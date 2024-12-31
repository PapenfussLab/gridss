package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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
			int indelCigarElementOffset,
		  	boolean isInAssemblyAnchor) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, offsetRemoteStart, offsetRemoteEnd, 0, 0, isInAssemblyAnchor);
		this.indel = indel;
		this.indelCigarElementOffset = indelCigarElementOffset;
	}
	public static IndelEvidence create(SAMEvidenceSource source, SAMRecord record, int indelCigarElementOffset) {
		boolean isInAssemblyAnchor = isEntirelyContainedInAssemblyAnchor(record, new ChimericAlignment(record), null);
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
				preStartOffset, preEndOffset, unmappedStartOffset, unmappedEndOffset, postStartOffset, postEndOffset, indel, indelCigarElementOffset, isInAssemblyAnchor);
		IndelEvidence right = new IndelEvidence(source, record, location.remoteBreakpoint(),
				postStartOffset, postEndOffset, unmappedStartOffset, unmappedEndOffset, preStartOffset, preEndOffset, indel, indelCigarElementOffset, isInAssemblyAnchor);
		left.remote = right;
		right.remote = left;
		return left;
	}
	public static List<IndelEvidence> create(SAMEvidenceSource source, int minIndelSize, SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) return Collections.emptyList();
		if (CigarUtil.widthOfImprecision(record.getCigar()) > 0) {
			// not a real indel: this is a placeholder CIGAR for an unanchored breakend assembly
			return Collections.emptyList();
		}
		List<IndelEvidence> list = new ArrayList<IndelEvidence>(4);
		List<CigarElement> cl = record.getCigar().getCigarElements();
		// find all indels (excluding those at start/end)
		int indelStartOffset = 1;
		int indelSize = 0;
		while (indelStartOffset < cl.size()) {
			if (cl.get(indelStartOffset).getOperator().isIndelOrSkippedRegion() &&
					!cl.get(indelStartOffset - 1).getOperator().isIndelOrSkippedRegion()) {
				int indelEndOffset = indelStartOffset;
				while (indelEndOffset < cl.size() && cl.get(indelEndOffset).getOperator().isIndelOrSkippedRegion()) {
					indelSize += cl.get(indelEndOffset).getLength();
					indelEndOffset++;
				}
				if (indelSize >= minIndelSize && indelEndOffset < cl.size()) {
					IndelEvidence left = create(source, record, indelStartOffset);
					IndelEvidence right = left.asRemote();
					list.add(left);
					list.add(right);
				}
				indelStartOffset = indelEndOffset;
				indelSize = 0;
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
			return scoreAssembly();
		}
		CigarElement e = indel.get(0);
		for (int i = 1; i < indel.size(); i++) {
			if (e.getLength() < indel.get(i).getLength()) {
				e = indel.get(i);
			}
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreIndel(getEvidenceSource().getMetrics(), this, e.getOperator(), e.getLength(), getLocalMapq());
	}
	@Override
	public float getBreakpointQual() {
		return getBreakendQual();
	}
	private float scoreAssembly() {
		AssemblyAttributes attr = new AssemblyAttributes(getSAMRecord());
		int pos = getBreakendAssemblyContigOffset();
		int rp = attr.getSupportingReadCount(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null, source.getContext());
		double rpq = attr.getSupportingQualScore(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null, source.getContext());
		int sc = attr.getSupportingReadCount(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null, source.getContext());
		double scq = attr.getSupportingQualScore(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null, source.getContext());
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreAssembly(this,
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
		return source.getContext().getEvidenceIDGenerator().getEvidenceID(this);
	}
	/**
	 * Identifies which indel in the read this evidence corresponds to.
	 * @return zero-based offset in the read CIGAR operator list of this indel
	 */
	public int getIndelCigarOffset() {
		return indelCigarElementOffset;
	}
	@Override
	public boolean isReference() {
		// TODO Auto-generated method stub
		return false;
	}
}
