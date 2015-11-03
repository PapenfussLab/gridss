package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;

import com.google.common.collect.ImmutableList;

/**
 * Assembly spanning a small indel
 */
public class SmallIndelSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence {
	private final SAMRecord assembly;
	/**
	 * Creates the forward breakpoint of the event spanned by the assembly
	 * @param source
	 * @param assembly
	 */
	public SmallIndelSAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly) {
		super(source, createFwdAnchored(assembly), ImmutableList.of(createFwdRealigned(assembly)));
		this.assembly = assembly;
		getRemoteSAMRecord().setMappingQuality(getSAMRecord().getMappingQuality());
		assert(sanityCheck());
	}
	/**
	 * Creates the backward breakpoint of the event spanned by the assembly
	 * @param assembly assembly
	 */
	private SmallIndelSAMRecordAssemblyEvidence(SmallIndelSAMRecordAssemblyEvidence assembly) {
		super(assembly.getEvidenceSource(), createBwdAnchored(assembly.getBackingRecord()), ImmutableList.of(createBwdRealigned(assembly.getBackingRecord())));
		this.assembly = assembly.assembly;
		getRemoteSAMRecord().setMappingQuality(getSAMRecord().getMappingQuality());
	}
	/**
	 * Emulates the anchored portion as if the assembly was a split read mapping
	 * @param indelAssembly assembly spanning indel breakpoint
	 * @param 
	 * @return anchored breakend alignment
	 */
	private static SAMRecord createFwdAnchored(final SAMRecord r) {
		SAMRecord read = createFwdSoftClipped(r);
		read.setReadName(getEvidenceID(r, BreakendDirection.Forward));
		return read;
	}
	private static SAMRecord createBwdAnchored(final SAMRecord r) {
		SAMRecord read = createBwdSoftClipped(r);
		read.setReadName(getEvidenceID(r, BreakendDirection.Backward));
		return read;
	}
	/**
	 * Emulates a split read realignment as if the assembly was a split read mapping
	 * @param indelAssembly assembly spanning indel breakpoint
	 * @return breakend split read realignment
	 */
	private static SAMRecord createFwdRealigned(final SAMRecord r) {
		SAMRecord read = createBwdSoftClipped(r);
		SAMRecordUtil.trimSoftClips(read, CigarUtil.readLength(CigarUtil.splitAtLargestIndel(r.getCigar().getCigarElements()).get(0)), 0);
		//read.setReadName(BreakpointFastqEncoding.getEncodedFastqID(r.getReferenceIndex(), createFwdAnchored(r).getAlignmentStart(), 0, getEvidenceID(r, BreakendDirection.Forward)));
		read.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY); // hack to force remote portion of assembly to not be mapq filtered
		return read;
	}
	private static SAMRecord createBwdRealigned(final SAMRecord r) {
		SAMRecord read = createFwdSoftClipped(r);
		// Realigned read should not include bases mapped to the anchored read
		SAMRecordUtil.trimSoftClips(read, 0, CigarUtil.readLength(CigarUtil.splitAtLargestIndel(r.getCigar().getCigarElements()).get(2)));
		//read.setReadName(BreakpointFastqEncoding.getEncodedFastqID(r.getReferenceIndex(), createBwdAnchored(r).getAlignmentStart(), 0, getEvidenceID(r, BreakendDirection.Forward)));
		read.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY); // hack to force remote portion of assembly to not be mapq filtered
		return read;
	}
	/**
	 * Creates a softclip mapping of the start of the read
	 * @param indelAssembly assembly spanning indel breakpoint
	 * @param 
	 * @return anchored breakend alignment
	 */
	private static SAMRecord createFwdSoftClipped(final SAMRecord r) {
		List<CigarElement> cigar = CigarUtil.splitAtLargestIndel(r.getCigar().getCigarElements()).get(0);
		int clipLength = r.getReadLength() - CigarUtil.readLength(cigar);
		if (clipLength > 0) {
			cigar.add(new CigarElement(clipLength, CigarOperator.SOFT_CLIP));
		}
		SAMRecord read = SAMRecordUtil.clone(r);
		read.setCigar(new Cigar(cigar));
		read.setAttribute(SamTags.ASSEMBLY_DIRECTION, BreakendDirection.Forward.toChar());
		read.setAttribute(SamTags.SPANNING_ASSEMBLY, (char)'y');
		return read;
	}
	/**
	 * Emulates a split read realignment as if the assembly was a split read mapping
	 * @param indelAssembly assembly spanning indel breakpoint
	 * @return breakend split read realignment
	 */
	private static SAMRecord createBwdSoftClipped(final SAMRecord r) {
		List<List<CigarElement>> split = CigarUtil.splitAtLargestIndel(r.getCigar().getCigarElements());
		List<CigarElement> cigar = split.get(2);
		int clipLength = r.getReadLength() - CigarUtil.readLength(cigar);
		if (clipLength > 0) {
			cigar.add(0, new CigarElement(clipLength, CigarOperator.SOFT_CLIP));
		}
		int bwdStartPosition = r.getAlignmentStart() + CigarUtil.referenceLength(split.get(0)) + CigarUtil.referenceLength(split.get(1));
		SAMRecord read = SAMRecordUtil.clone(r);
		read.setAlignmentStart(bwdStartPosition);
		read.setCigar(new Cigar(cigar));
		read.setAttribute(SamTags.ASSEMBLY_DIRECTION, BreakendDirection.Backward.toChar());
		read.setAttribute(SamTags.SPANNING_ASSEMBLY, (char)'y');
		return read;
	}
	@Override
	public SAMRecord getBackingRecord() {
		return assembly;
	}
	@Override
	public SmallIndelSAMRecordAssemblyEvidence asRemote() {
		if (getBreakendSummary().direction == BreakendDirection.Forward) {
			return new SmallIndelSAMRecordAssemblyEvidence(this);
		} else {
			return new SmallIndelSAMRecordAssemblyEvidence(this.getEvidenceSource(), this.assembly);
		}
	}
	@Override
	public String getEvidenceID() {
		return getEvidenceID(assembly, getBreakendSummary().direction);
	}
	private static String getEvidenceID(SAMRecord r, BreakendDirection direction) {
		String prefix = direction == BreakendDirection.Forward ? "f" : "b";
		return prefix + r.getReadName();
	}
	@Override
	public boolean isSpanningAssembly() {
		return true;
	}
	@Override
	public SmallIndelSAMRecordAssemblyEvidence annotateAssembly() {
		super.annotateAssembly();
		getRemoteSAMRecord().setMappingQuality(getLocalMapq());
		return this;
	}
	private boolean sanityCheck() {
		List<List<CigarElement>> split = CigarUtil.splitAtLargestIndel(assembly.getCigar().getCigarElements());
		if (getBreakendSummary()  != null) {
			assert(getBreakendSummary() instanceof BreakpointSummary);
			assert(getBreakendSummary().direction == BreakendDirection.Forward);
			assert(getBreakendSummary().direction2 == BreakendDirection.Backward);
			assert(split.get(0).size() > 0);
			assert(split.get(1).size() > 0);
			assert(split.get(2).size() > 0);
		} else {
			assert(isReferenceAssembly());
		}
		return true;
	}
}
