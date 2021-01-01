package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.NotImplementedException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Chimeric alignment based support for a structural variant
 * @author Daniel Cameron
 *
 */
public class SplitReadEvidence extends SingleReadEvidence implements DirectedBreakpoint {
	private static final Log log = Log.getInstance(SplitReadEvidence.class);
	private ChimericAlignment remoteAlignment;
	protected SplitReadEvidence(SAMEvidenceSource source, SAMRecord record, BreakendSummary location,
			int offsetLocalStart, int offsetLocalEnd,
			int offsetUnmappedStart, int offsetUnmappedEnd,
			int offsetRemoteStart, int offsetRemoteEnd,
			ChimericAlignment remoteAlignment,
			int localUnanchoredWidth, int remoteUnanchoredWidth,
			boolean isInAssemblyAnchor) {
		super(source, record, location, offsetLocalStart, offsetLocalEnd, offsetUnmappedStart, offsetUnmappedEnd, offsetRemoteStart, offsetRemoteEnd, localUnanchoredWidth, remoteUnanchoredWidth, isInAssemblyAnchor);
		this.remoteAlignment = remoteAlignment;
	}
	public static List<SplitReadEvidence> create(SAMEvidenceSource source, SAMRecord record) {
		if (record.getReadUnmappedFlag() || record.getCigar() == null) return Collections.emptyList();
		List<ChimericAlignment> aln = ChimericAlignment.getChimericAlignments(record);
		if (aln.isEmpty()) return Collections.emptyList();
		if (record.getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP
				|| record.getCigar().getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
			record = record.deepCopy();
			SAMRecordUtil.hardClipToN(record);
		}
		List<SplitReadEvidence> list = new ArrayList<>(2);
		ChimericAlignment chim = new ChimericAlignment(record);
		int offset = SAMRecordUtil.getFirstAlignedBaseReadOffset(record);
		ChimericAlignment pre = aln.stream().filter(ca -> ca.getFirstAlignedBaseReadOffset() < offset).max(ChimericAlignment.ByReadOffset).orElse(null);
		ChimericAlignment post = aln.stream().filter(ca -> ca.getFirstAlignedBaseReadOffset() > offset).min(ChimericAlignment.ByReadOffset).orElse(null);
		SAMSequenceDictionary dict = source != null ? source.getContext().getDictionary() : record.getHeader().getSequenceDictionary();
		// Read is AAAXBBBYCCC
		// we are alignment B
		// pre is A
		// post is C
		// X,Y are unaligned bases
		final int rl = record.getReadLength();
		int startOffset = chim.getFirstAlignedBaseReadOffset();
		int endOffset = chim.getLastAlignedBaseReadOffset() + 1;
		if (pre != null) {
			int overlap = pre.getLastAlignedBaseReadOffset() + 1 - startOffset;
			BreakpointSummary bs = new BreakpointSummary(withOverlap(chim.predecessorBreakend(dict), overlap), withOverlap(pre.successorBreakend(dict), overlap));
			int preStartOffset = 0; // ignore the actual alignment and just go out to the end of the read so we can assemble across multiple breakpoints
			int preEndOffset = pre.getLastAlignedBaseReadOffset() + 1;
			if (record.getReadNegativeStrandFlag()) {
				list.add(new SplitReadEvidence(source, record, bs,
						rl - (INCLUDE_CLIPPED_ANCHORING_BASES ? record.getReadLength() : endOffset), rl - startOffset,
						rl - startOffset, rl - preEndOffset,
						rl - preEndOffset, rl - preStartOffset,
						pre,
						CigarUtil.widthOfImprecision(chim.cigar), CigarUtil.widthOfImprecision(pre.cigar),
						isEntirelyContainedInAssemblyAnchor(record, chim, pre)));
			} else {
				list.add(new SplitReadEvidence(source, record, bs,
						startOffset, INCLUDE_CLIPPED_ANCHORING_BASES ? record.getReadLength() : endOffset,
						preEndOffset, startOffset,
						preStartOffset, preEndOffset,
						pre,
						CigarUtil.widthOfImprecision(chim.cigar), CigarUtil.widthOfImprecision(pre.cigar),
						isEntirelyContainedInAssemblyAnchor(record, chim, pre)));
			}
		}
		if (post != null) {
			int overlap = endOffset - post.getFirstAlignedBaseReadOffset();
			BreakpointSummary bs = new BreakpointSummary(withOverlap(chim.successorBreakend(dict), overlap), withOverlap(post.predecessorBreakend(dict), overlap));
			int postStartOffset = post.getFirstAlignedBaseReadOffset();
			int postEndOffset = rl;
			if (record.getReadNegativeStrandFlag()) {
				list.add(new SplitReadEvidence(source, record, bs,
						rl - endOffset, rl - (INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : startOffset),
						rl - postStartOffset, rl - endOffset, 
						rl - postEndOffset, rl - postStartOffset, 
						post,
						CigarUtil.widthOfImprecision(chim.cigar), CigarUtil.widthOfImprecision(post.cigar),
						isEntirelyContainedInAssemblyAnchor(record, chim, post)));
			} else {
				list.add(new SplitReadEvidence(source, record, bs,
					INCLUDE_CLIPPED_ANCHORING_BASES ? 0 : startOffset, endOffset,
					endOffset, postStartOffset,
					postStartOffset, postEndOffset,
					post,
					CigarUtil.widthOfImprecision(chim.cigar), CigarUtil.widthOfImprecision(post.cigar),
					isEntirelyContainedInAssemblyAnchor(record, chim, post)));
			}
		}
		return list;
	}
	/**
	 * Adjusts the breakend bounds assuming an overalignment due to microhomology 
	 * @param bs nominal breakend position
	 * @param overlap number of bases overlapping with the adjacent chimeric alignment
	 * @return Breakend position factoring in microhomology
	 */
	private static BreakendSummary withOverlap(BreakendSummary bs, int overlap) {
		if (overlap <= 0) return bs;
		return new BreakendSummary(bs.referenceIndex, bs.direction, bs.nominal,
				bs.start - (bs.direction == BreakendDirection.Forward ? overlap : 0),
				bs.end + (bs.direction == BreakendDirection.Backward ? overlap : 0));
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return (BreakpointSummary)super.getBreakendSummary();
	}

	@Override
	public int getRemoteMapq() {
		return remoteAlignment.mapq;
	}
	@Override
	public float getBreakendQual() {
		return getBreakpointQual();
	}
	@Override
	public float getBreakpointQual() {
		if (AssemblyAttributes.isAssembly(getSAMRecord())) {
			return scoreAssembly();
		}
		int softClipLength = getBreakendSequence().length;
		if (getSAMRecord().getSupplementaryAlignmentFlag()) {
			ChimericAlignment caThis = new ChimericAlignment(getSAMRecord());
			// The first record should be the primary
			ChimericAlignment caPrimary = ChimericAlignment.getChimericAlignments(getSAMRecord().getStringAttribute(SAMTag.SA.name())).get(0);
			ChimericAlignment osaPrimary = caPrimary;
			String osa = (String)getSAMRecord().getTransientAttribute("OSA");
			if (osa != null) {
				// Use the original SA if it exists - the primary might be considered unmapped and not in the SA tag
				osaPrimary = ChimericAlignment.getChimericAlignments(osa).get(0);
			}
			if (!caPrimary.equals(osaPrimary)) {
				// We've stripped the primary alignment.
				// Since the primary could be in the middle, we'll just take the shorter of the primary soft clips
				int myOffset = getLocalChimericAlignmentReadOffset();
				int primaryOffset = osaPrimary.getFirstAlignedBaseReadOffset();
				int remoteOffset = getRemoteChimericAlignmentReadOffset();
				if (Math.min(myOffset, remoteOffset) <= primaryOffset && Math.max(myOffset, remoteOffset) >= primaryOffset) {
					return Math.min(SAMRecordUtil.getStartClipLength(osaPrimary.cigar.getCigarElements()), SAMRecordUtil.getEndClipLength(osaPrimary.cigar.getCigarElements()));
				}
				caPrimary = osaPrimary;
			}
			BreakendDirection primaryDirectionTowardThis = caThis.getFirstAlignedBaseReadOffset() < caPrimary.getFirstAlignedBaseReadOffset() ^ caPrimary.isNegativeStrand ? BreakendDirection.Backward : BreakendDirection.Forward;
			softClipLength = SAMRecordUtil.getSoftClipLength(caPrimary.cigar.getCigarElements(), primaryDirectionTowardThis);
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreSplitRead(getEvidenceSource().getMetrics(),
				softClipLength,
				getLocalMapq(), getRemoteMapq());
	}
	private float scoreAssembly() {
		AssemblyAttributes attr = new AssemblyAttributes(getSAMRecord());
		int pos = getBreakendAssemblyContigOffset();
		int rp = attr.getSupportingReadCount(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null);
		double rpq = attr.getSupportingQualScore(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair), null);
		int sc = attr.getSupportingReadCount(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null);
		double scq = attr.getSupportingQualScore(pos, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null);
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreAssembly(
				rp, rpq,
				sc, scq,
				getLocalMapq(),
				getRemoteMapq());
	}
	@Override
	public DirectedBreakpoint asRemote() {
		throw new NotImplementedException("asRemote() should no longer be required");
	}
	@Override
	protected String getUncachedEvidenceID() {
		return source.getContext().getEvidenceIDGenerator().getEvidenceID(this);
	}
	@Override
	public String getRemoteEvidenceID() {
		SAMRecord remote = this.getSAMRecord().deepCopy();
		remote.setReferenceName(remoteAlignment.rname);
		remote.setAlignmentStart(remoteAlignment.pos);
		remote.setReadUnmappedFlag(false);
		remote.setReadNegativeStrandFlag(remoteAlignment.isNegativeStrand);
		remote.setCigar(remoteAlignment.cigar);
		remote.setAttribute(SAMTag.SA.name(), new ChimericAlignment(this.getSAMRecord()).toString());
		SplitReadEvidence remoteEvidence = SplitReadEvidence.create(source, remote).get(0);
		return source.getContext().getEvidenceIDGenerator().getEvidenceID(remoteEvidence);
	}
	@Override
	public boolean involvesPrimaryReadAlignment() {
		return super.involvesPrimaryReadAlignment()
			// the first record in the SA tag should be the primary read alignment 
			|| ChimericAlignment.getChimericAlignments(getSAMRecord()).get(0).equals(remoteAlignment);
	}

	/**
	 * Evidence provides support for no structural variant call
	 * This can be due to sequence homology allowing the alignment
	 * at either breakend to be placed on the opposite breakend
	 */
	@Override
	public boolean isReference() {
		if (!isBreakendExact()) return false;
		if (getUntemplatedSequence().length() > 0) return false;
		BreakendSummary location = getBreakendSummary();
		int anchorHomLen = location.nominal - location.start;
		int remoteHomLen = location.end - location.nominal;
		if (location.direction == BreakendDirection.Backward) {
			int tmp = anchorHomLen;
			anchorHomLen = remoteHomLen;
			remoteHomLen = tmp;
		}
		ChimericAlignment localAlignment = new ChimericAlignment(getSAMRecord()); 
		int localBases = localAlignment.getLastAlignedBaseReadOffset() - localAlignment.getFirstAlignedBaseReadOffset() + 1;
		int remoteBases = remoteAlignment.getLastAlignedBaseReadOffset() - remoteAlignment.getFirstAlignedBaseReadOffset() + 1;
		return anchorHomLen >= localBases || remoteHomLen >= remoteBases;
	}
	/**
	 * Offset of the read alignment
	 * @return
	 */
	public int getRemoteChimericAlignmentReadOffset() {
		return remoteAlignment.getFirstAlignedBaseReadOffset();
	}

}
