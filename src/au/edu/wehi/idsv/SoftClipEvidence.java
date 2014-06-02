package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.SequenceUtil;

import java.nio.charset.StandardCharsets;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.VcfAttributes;

public class SoftClipEvidence implements DirectedBreakpoint {
	private final ProcessingContext processContext;
	private final SAMRecord record;
	private final BreakendSummary location;
	private final RealignedBreakpoint realigned;
	public SoftClipEvidence(ProcessingContext processContext, BreakendDirection direction, SAMRecord record, SAMRecord realigned) {
		if (record == null) throw new IllegalArgumentException("record is null");
		if (direction == null) throw new IllegalArgumentException("direction is null");
		if (record.getReadUnmappedFlag()) throw new IllegalArgumentException(String.format("record %s is unmapped", record.getReadName()));
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE ) throw new IllegalArgumentException(String.format("record %s missing sequence information", record.getReadName()));
		this.processContext = processContext;
		this.record = record;
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart();
		BreakendSummary local = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos, new EvidenceMetrics());		
		local.evidence.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 1);
		//local.evidence.set(EvidenceAttributes.SOFT_CLIP_MAX_LENGTH, getSoftClipLength());
		local.evidence.set(VcfAttributes.SOFT_CLIP_TOTAL_LENGTH, getSoftClipLength(direction, record));
		if (realigned != null && !realigned.getReadUnmappedFlag()) {
			this.realigned = new RealignedBreakpoint(local, realigned);
			this.location = this.realigned.getBreakpointSummary();
		} else {
			this.location = local;
			this.realigned = null;
		}
		if (getSoftClipLength() == 0) {
			throw new IllegalArgumentException(String.format("record %s is not %s soft clipped", record.getReadName(), direction));
		}
	}
	public static int getSoftClipLength(BreakendDirection direction, SAMRecord record) {
		return direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(record) : SAMRecordUtil.getStartSoftClipLength(record); 
	}
	public SoftClipEvidence(ProcessingContext processContext, BreakendDirection direction, SAMRecord record) {
		this(processContext, direction, record, null);
	}
	public SoftClipEvidence(SoftClipEvidence evidence, SAMRecord realigned) {
		this(evidence.processContext, evidence.location.direction, evidence.record, realigned);
	}
	@Override
	public String getEvidenceID() {
		// need read name, breakpoint direction & which read in pair
		String readNumber = record.getReadPairedFlag() ? record.getFirstOfPairFlag() ? "/1" : "/2" : "";
		return String.format("%s%s%s", location.direction == BreakendDirection.Forward ? "f" : "b", record.getReadName(), readNumber);
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	@Override
	public byte[] getBreakpointSequence() {
		return location.direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBases(record) : SAMRecordUtil.getStartSoftClipBases(record);
	}
	@Override
	public byte[] getBreakpointQuality() {
		return location.direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBaseQualities(record) : SAMRecordUtil.getStartSoftClipBaseQualities(record);
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public int getSoftClipLength() {
		return getSoftClipLength(location.direction, record); 
	}
	
	/**
	 * Number of unmapped bases at the breakpoint 
	 * @return Number of unmapped bases at the breakpoint
	 */
	public int getUntemplatedSequenceLength() {
		if (location instanceof BreakpointSummary) {
			return realigned.getInsertedSequenceLength();
		} else {
			return getSoftClipLength();
		}
	}
	/**
	 * Sequence of untemplated bases
	 * @return
	 */
	public String getUntemplatedSequence() {
		if (location instanceof BreakpointSummary) {
			return realigned.getInsertedSequence();
		} else {
			return new String(getBreakpointSequence(), StandardCharsets.US_ASCII);
		}
	}
	/**
	 * 0-100 scaled percentage identity of mapped read bases.
	 * @return percentage of reference-aligned bases that match reference 
	 */
	public float getAlignedPercentIdentity() {
		// final byte[] referenceBases = refSeq.get(sequenceDictionary.getSequenceIndex(rec.getReferenceName())).getBases();
        // rec.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(rec, referenceBases, 0, bisulfiteSequence));
        //if (rec.getBaseQualities() != SAMRecord.NULL_QUALS) {
        // rec.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(rec, referenceBases, 0, bisulfiteSequence));
        Integer nm = record.getIntegerAttribute(SAMTag.NM.name());
		if (nm != null) {
			int refBasesToConsider = record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record) - SAMRecordUtil.getEndSoftClipLength(record); 
			int refBaseMatches = refBasesToConsider - nm + SequenceUtil.countInsertedBases(record) + SequenceUtil.countDeletedBases(record); 
			return 100.0f * refBaseMatches / (float)refBasesToConsider;
		}
		String md = record.getStringAttribute(SAMTag.MD.name());
		if (StringUtils.isNotEmpty(md)) {
			// Socrates handles this: should we? Which aligners write MD but not NM?
			throw new RuntimeException("Sanity Check Failure: Not Yet Implemented: calculation from reads with MD tag but not NM tag as per Socrates implementation");
		}
		throw new IllegalStateException(String.format("Read %s missing NM tag", record.getReadName()));
	}
	public float getAverageClipQuality() {
		float total = 0;
		byte[] qual = getBreakpointQuality();
		if (qual == null) return 0;
		for (int i = 0; i < qual.length; i++) {
			total += qual[i]; 
		}
		return total / qual.length;
	}
	public int getMappingQuality() {
		return record.getMappingQuality();
	}
}
