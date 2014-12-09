package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.SequenceUtil;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class SoftClipEvidence implements DirectedEvidence {
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final SAMRecord record;
	private final BreakendSummary location;
	/**
	 * Lazily computed evidence identifier
	 */
	private String evidenceID = null;
	public static SoftClipEvidence create(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		return create(processContext, source, direction, record, null);
	}
	public static SoftClipEvidence create(SoftClipEvidence evidence, SAMRecord realigned) {
		return create(evidence.processContext, evidence.source, evidence.location.direction, evidence.record, realigned);
	}
	public static SoftClipEvidence create(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record, SAMRecord realigned) {
		if (record == null) throw new IllegalArgumentException("record is null");
		if (direction == null) throw new IllegalArgumentException("direction is null");
		if (record.getReadUnmappedFlag()) throw new IllegalArgumentException(String.format("record %s is unmapped", record.getReadName()));
		if (record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE ) throw new IllegalArgumentException(String.format("record %s missing sequence information", record.getReadName()));
		SoftClipEvidence result = null;
		if (realigned != null && !realigned.getReadUnmappedFlag() && processContext.getRealignmentParameters().realignmentPositionUnique(realigned)) {
			try {
				result = new RealignedSoftClipEvidence(processContext, source, direction, record, realigned);
			} catch (CloneNotSupportedException e) {
				throw new RuntimeException(e);
			}
		} else {
			result = new SoftClipEvidence(processContext, source, direction, record);
		}
		if (result.getSoftClipLength() == 0) {
			throw new IllegalArgumentException(String.format("record %s is not %s soft clipped", record.getReadName(), direction));
		}
		return result;
	}
	protected SoftClipEvidence(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record) {
		this.processContext = processContext;
		this.source = source;
		this.record = record;
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart();
		this.location = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos);
	}
	public static int getSoftClipLength(BreakendDirection direction, SAMRecord record) {
		return direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(record) : SAMRecordUtil.getStartSoftClipLength(record); 
	}
	public static String getEvidenceID(BreakendDirection direction, SAMRecord softClip) {
		// need read name, breakpoint direction & which read in pair
		String readNumber = softClip.getReadPairedFlag() ? softClip.getFirstOfPairFlag() ? "/1" : "/2" : "";
		return String.format("%s%s%s", direction == BreakendDirection.Forward ? "f" : "b", softClip.getReadName(), readNumber);
	}
	@Override
	public String getEvidenceID() {
		if (evidenceID == null) {
			evidenceID = getEvidenceID(location.direction, record);
		}
		return evidenceID;
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return location;
	}
	@Override
	public byte[] getBreakendSequence() {
		return location.direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBases(record) : SAMRecordUtil.getStartSoftClipBases(record);
	}
	@Override
	public byte[] getBreakendQuality() {
		return location.direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBaseQualities(record) : SAMRecordUtil.getStartSoftClipBaseQualities(record);
	}
	public SAMRecord getSAMRecord() {
		return this.record;
	}
	public int getSoftClipLength() {
		return getSoftClipLength(location.direction, record); 
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
		byte[] qual = getBreakendQuality();
		if (qual == null) return 0;
		for (int i = 0; i < qual.length; i++) {
			total += qual[i]; 
		}
		return total / qual.length;
	}
	public int getMappingQuality() {
		return record.getMappingQuality();
	}
	public SAMEvidenceSource getEvidenceSource() {
		return source;
	}
	@Override
	public int getLocalMapq() {
		return record.getMappingQuality();
	}
	@Override
	public int getLocalBaseLength() {
		return record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record) - SAMRecordUtil.getEndSoftClipLength(record);
	}
	@Override
	public int getLocalMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(record);
	}
	@Override
	public int getLocalTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(record);
	}
	@Override
	public String toString() {
		return "SoftClip len=" + getSoftClipLength() + " " + getBreakendSummary().toString() + " " + getSAMRecord().getReadName();
	}
	/**
	 * Determines whether this evidence provides support for a putative SV
	 * @param p soft clip parameters
	 * @param rm metrics
	 * @return true if the soft clip provides support, false otherwise
	 */
	public boolean meetsEvidenceCritera(SoftClipParameters p) {
		return getMappingQuality() >= p.minReadMapq
				&& getSoftClipLength() >= p.minLength
				&& getAlignedPercentIdentity() >= p.minAnchorIdentity
				&& !isDovetailing()
				&& !p.adapters.isAdapterSoftClip(this);
	}
	/**
	 * Dovetailing reads do not support SVs, they are caused by fragment size
	 * less than read length
	 * 
	 * =======> <=======
	 * 
	 * 
	 * @param expectedOrientation
	 *            read pair orientation
	 * @return true if the soft clip is due to a fragment size smaller than the
	 *         read length
	 */
	public boolean isDovetailing() {
		if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()
				|| record.getReadUnmappedFlag())
			return false;
		return record.getMateReferenceIndex() == record.getReferenceIndex()
				&& Math.abs(record.getAlignmentStart()
						- record.getMateAlignmentStart()) <= Defaults.READ_PAIR_DOVETAIL_MARGIN
				// dovetails happen on the 3' end of the read for FR
				&& ((location.direction == BreakendDirection.Forward && !record
						.getReadNegativeStrandFlag()) || (location.direction == BreakendDirection.Backward && record
						.getReadNegativeStrandFlag()));
	}
	@Override
	public boolean isBreakendExact() {
		return true;
	}
}
