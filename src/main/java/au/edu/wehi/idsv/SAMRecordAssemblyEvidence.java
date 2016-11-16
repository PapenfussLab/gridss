package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.primitives.Bytes;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.sam.SplitIndel;
import au.edu.wehi.idsv.util.MathUtil;
import au.edu.wehi.idsv.vcf.VcfFilter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

public class SAMRecordAssemblyEvidence implements AssemblyEvidence {
	private static final Log log = Log.getInstance(SAMRecordAssemblyEvidence.class);
	public static final String COMPONENT_EVIDENCEID_SEPARATOR = " ";
	private final SAMRecord record;
	private final SAMRecord realignment;
	private final CompoundBreakendAlignment realignments;
	private final AssemblyEvidenceSource source;
	private final BreakendSummary breakend;
	private final boolean isExact;
	private List<SpanningSAMRecordAssemblyEvidence> indels = null;
	private Collection<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
	private HashSet<String> evidenceIDs = null;
	public SAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecordAssemblyEvidence assembly, List<SAMRecord> realigned) {
		this(source, assembly.getSAMRecord(), realigned);
		this.evidence = assembly.evidence;
		this.evidenceIDs = assembly.evidenceIDs;
	}
	public SAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly, List<SAMRecord> realigned) {
		this(source, assembly, realigned,
				calculateBreakend(assembly),
				calculateIsBreakendExact(assembly.getCigar()));
	}
	private SAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly, List<SAMRecord> realigned, BreakendSummary bs, boolean isExact) {
		this.source = source;
		this.record = assembly;
		this.isExact = isExact;
		this.realignments = new CompoundBreakendAlignment(source.getContext(), assembly, realigned);
		this.breakend = bs;
		this.realignment = realignments.getSimpleBreakendRealignment();
		SAMRecordUtil.pairReads(this.record, this.realignment);
	}
	public List<SAMRecordAssemblyEvidence> getSubsequentRealignments() {
		if (realignments.getBreakpointCount() <= 1) {
			return ImmutableList.of();
		}
		List<SAMRecordAssemblyEvidence> list = new ArrayList<SAMRecordAssemblyEvidence>(realignments.getBreakpointCount());
		for (Pair<SAMRecord, SAMRecord> pair : realignments.getSubsequentBreakpointAlignmentPairs()) {
			SAMRecord anchor = pair.getLeft();
			SAMRecord realign = pair.getRight();
			SAMRecordAssemblyEvidence be = AssemblyFactory.hydrate(getEvidenceSource(), anchor);
			SAMRecordAssemblyEvidence bp = AssemblyFactory.incorporateRealignment(getEvidenceSource().getContext(), be, ImmutableList.of(realign));
			list.add(bp);
		}
		return list;
	}
	/**
	 * Gets all indels spanned by this event
	 * @return
	 */
	public List<SpanningSAMRecordAssemblyEvidence> getSpannedIndels() {
		if (indels == null) {
			List<SplitIndel> splitindels = SplitIndel.getIndelsAsSplitReads(getSAMRecord());
			indels = new ArrayList<SpanningSAMRecordAssemblyEvidence>(splitindels.size() * 2);
			int i = 0;
			for (SplitIndel indel : splitindels) {
				indel = indel.onPositive();
				indel.leftAnchored.setAttribute(SamTags.ASSEMBLY_DIRECTION, BreakendDirection.Forward.toChar());
				indel.rightAnchored.setAttribute(SamTags.ASSEMBLY_DIRECTION, BreakendDirection.Backward.toChar());
				SpanningSAMRecordAssemblyEvidence left = new SpanningSAMRecordAssemblyEvidence(this, indel, i++);
				indels.add(left);
				indels.add(left.asRemote());
			}
		}
		return indels;
	}
	/**
	 * Determines whether the given record is part of the given assembly
	 *
	 * This method is a probabilistic method and it is possible for the record to return true
	 * when the record does not form part of the assembly breakend
	 *
	 * @param evidence
	 * @return true if the record is likely part of the breakend, false if definitely not
	 */
	public boolean isPartOfAssembly(DirectedEvidence e) {
		return getEvidenceIDs().contains(e.getEvidenceID());
	}
	/**
	 * Hydrates the given evidence back into the assembly evidence set
	 * @param e
	 */
	public SAMRecordAssemblyEvidence hydrateEvidenceSet(DirectedEvidence e) {
		if (isPartOfAssembly(e)) {
			evidence.add(e);
		}
		return this;
	}
	/**
	 * Hydrates the given evidence back into the assembly evidence set
	 * @param e
	 */
	public SAMRecordAssemblyEvidence hydrateEvidenceSet(Collection<? extends DirectedEvidence> evidence) {
		for (DirectedEvidence e : evidence) {
			hydrateEvidenceSet(e);
		}
		return this;
	}
	/**
	 * Annotates an assembly with a fully hydrated evidence set
	 */
	public SAMRecordAssemblyEvidence annotateAssembly() {
		int n = source.getContext().getCategoryCount();
		BreakendSummary breakendWithMargin = source.getContext().getVariantCallingParameters().withMargin(breakend);
		float[] rpQual = new float[n];
		float[] scQual = new float[n];
		float[] rQual = new float[n];
		float[] nsrpQual = new float[n];
		float[] nsscQual = new float[n];
		int[] rpCount = new int[n];
		int[] rpMaxLen = new int[n];
		int[] scCount = new int[n];
		int[] scLenMax = new int[n];
		int[] scLenTotal = new int[n];
		int[] rCount = new int[n];
		int[] nsrpCount = new int[n];
		int[] nsscCount = new int[n];
		int maxLocalMapq = 0;
		for (DirectedEvidence e : getEvidence()) {
			maxLocalMapq = Math.max(maxLocalMapq, e.getLocalMapq());
			int offset = ((SAMEvidenceSource)e.getEvidenceSource()).getSourceCategory();
			float qual = e.getBreakendQual();
			if (e instanceof NonReferenceReadPair) {
				rpCount[offset]++;
				rpQual[offset] += qual;
				rpMaxLen[offset] = Math.max(rpMaxLen[offset], ((NonReferenceReadPair)e).getNonReferenceRead().getReadLength());
				if (breakendWithMargin != null && !breakendWithMargin.overlaps(e.getBreakendSummary())) {
					nsrpCount[offset]++;
					nsrpQual[offset] += qual;
				}
			}
			if (e instanceof SoftClipEvidence) {
				scCount[offset]++;
				scQual[offset] += qual;
				int clipLength = ((SoftClipEvidence)e).getSoftClipLength();
				scLenMax[offset] = Math.max(scLenMax[offset], clipLength);
				scLenTotal[offset] += clipLength;
				if (breakendWithMargin != null && !breakendWithMargin.overlaps(e.getBreakendSummary())) {
					nsscCount[offset]++;
					nsscQual[offset] += qual;
				}
				if (e instanceof RemoteEvidence) {
					rCount[offset]++;
					rQual[offset] += qual;
				}
			}
		}
		boolean untrackEvidence = !source.getContext().getConfig().getAssembly().trackEvidenceID;
		annotateSAMRecord(getSAMRecord(), rpQual, scQual, rQual, nsrpQual, nsscQual, rpCount, rpMaxLen, scCount, scLenMax, scLenTotal, rCount, nsrpCount, nsscCount, untrackEvidence);
		annotateSAMRecord(getBackingRecord(), rpQual, scQual, rQual, nsrpQual, nsscQual, rpCount, rpMaxLen, scCount, scLenMax, scLenTotal, rCount, nsrpCount, nsscCount, untrackEvidence);
		annotateSAMRecord(getRemoteSAMRecord(), rpQual, scQual, rQual, nsrpQual, nsscQual, rpCount, rpMaxLen, scCount, scLenMax, scLenTotal, rCount, nsrpCount, nsscCount, untrackEvidence);
		// TODO: proper mapq model
		getSAMRecord().setMappingQuality(maxLocalMapq);
		getBackingRecord().setMappingQuality(maxLocalMapq);
		if (evidenceIDs != null && evidence != null && evidenceIDs.size() != evidence.size()) {
			log.warn(String.format("Found only %d of %d reads supporting assembly %s", evidence.size(), evidenceIDs.size(), getEvidenceID()));
		}
		if (untrackEvidence) {
			evidenceIDs = null;
			evidence = null;
		}
		if (getLocalMapq() < source.getContext().getConfig().minMapq) {
			log.warn(String.format("Sanity check failure: %s has mapq below minimum", getEvidenceID()));
		}
		return this;
	}
	private static void annotateSAMRecord(SAMRecord r, float[] rpQual, float[] scQual,
			float[] rQual, float[] nsrpQual, float[] nsscQual, int[] rpCount, int[] rpMaxLen,
			int[] scCount, int[] scLenMax, int[] scLenTotal, int[] rCount,
			int[] nsrpCount, int[] nsscCount, boolean untrackEvidence) {
		r.setAttribute(SamTags.ASSEMBLY_READPAIR_COUNT, rpCount);
		r.setAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX, rpMaxLen);
		r.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT, scCount);
		r.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_REMOTE_COUNT, rCount);
		r.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT, nsrpCount);
		r.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT, nsscCount);
		r.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, scLenMax);
		r.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, scLenTotal);
		r.setAttribute(SamTags.ASSEMBLY_READPAIR_QUAL, rpQual);
		r.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL, scQual);
		r.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_REMOTE_QUAL, rQual);
		r.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL, nsrpQual);
		r.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL, nsscQual);
		if (untrackEvidence) {
			r.setAttribute(SamTags.EVIDENCEID, null);
		}
	}
	@Override
	public Collection<String> getEvidenceIDs() {
		if (evidenceIDs == null) {
			String encoded = getBackingRecord().getStringAttribute(SamTags.EVIDENCEID);
			if (encoded == null) {
				throw new IllegalStateException("Unable to get constituent evidenceIDs from assembly with evidence tracking disabled");
			}
			String[] ids = encoded.split(COMPONENT_EVIDENCEID_SEPARATOR);
			evidenceIDs = new HashSet<String>(Arrays.asList(ids));
			evidenceIDs.remove("");
		}
		return evidenceIDs;
	}
	public Collection<DirectedEvidence> getEvidence() {
		Collection<String> ids = getEvidenceIDs();
		if (evidence.size() < ids.size()) {
			if (this instanceof RealignedRemoteSAMRecordAssemblyEvidence) {
				// We don't expect rehydration of short SC, OEA, or alternatively mapped evidence at remote position 
			} else {
				log.debug(String.format("Expected %d, found %d support for %s %s%s", ids.size(), evidence.size(), getEvidenceID(), debugEvidenceMismatch(), debugEvidenceIDsMssingFromEvidence()));
			}
		}
		if (evidence.size() > ids.size()) {
			// Don't throw exception as the user can't actually do anything about this
			// Just continue with as correct results as we can manage
			log.debug(String.format("Expected %d, found %d support for %s %s%s", ids.size(), evidence.size(), getEvidenceID(), debugEvidenceMismatch(), debugEvidenceIDsMssingFromEvidence()));
			//log.warn("Hash collision has resulted in evidence being incorrectly attributed to assembly " + getEvidenceID());
		}
		return evidence;
	}
	private String debugEvidenceMismatch() {
		StringBuilder sb = new StringBuilder();
		for (String s : debugEvidenceIDsMssingFromEvidence()) {
			sb.append(" -");
			sb.append(s);
		}
		for (String s : debugEvidenceNotInAssembly()) {
			sb.append(" +");
			sb.append(s);
		}
		return sb.toString();
	}
	private List<String> debugEvidenceIDsMssingFromEvidence() {
		List<String> result = new ArrayList<String>();
		for (String s : getEvidenceIDs()) {
			boolean found = false;
			for (DirectedEvidence e : evidence) {
				if (e.getEvidenceID().equals(s)) {
					found = true;
					break;
				}
			}
			if (!found) {
				result.add(s);
			}
		}
		return result;
	}
	private List<String> debugEvidenceNotInAssembly() {
		List<String> result = new ArrayList<String>();
		for (DirectedEvidence e : evidence) {
			boolean found = false;
			for (String s : getEvidenceIDs()) {
				if (e.getEvidenceID().equals(s)) {
					found = true;
					break;
				}
			}
			if (!found) {
				result.add(e.getEvidenceID());
			}
		}
		return result;
	}
	/**
	 * Initial direction of assembly breakend.
	 * @return
	 */
	public BreakendDirection getBreakendDirection() {
		return getBreakendDirection(record);
	}
	private static BreakendDirection getBreakendDirection(SAMRecord record) {
		Character c = (Character)record.getAttribute(SamTags.ASSEMBLY_DIRECTION);
		if (c == null) return null;
		return BreakendDirection.fromChar((char)c);
	}
	private static BreakendSummary calculateBreakend(SAMRecord record) {
		BreakendDirection direction = getBreakendDirection(record);
		if (direction == null) return null;
		if (SAMRecordUtil.getSoftClipLength(record, direction) == 0) return null;
		int beStart;
		int beEnd;
		if (!calculateIsBreakendExact(record.getCigar())) {
			beStart = record.getAlignmentStart();
			beEnd = record.getAlignmentEnd();
		} else if (direction == BreakendDirection.Forward) {
			beStart = record.getAlignmentEnd();
			beEnd = record.getAlignmentEnd();
		} else {
			beStart = record.getAlignmentStart();
			beEnd = record.getAlignmentStart();
		}
		return new BreakendSummary(record.getReferenceIndex(), direction, MathUtil.average(beStart, beEnd), beStart, beEnd);
	}
	private static boolean calculateIsBreakendExact(Cigar cigar) {
		for (CigarElement ce : cigar.getCigarElements()) {
			if (ce.getOperator() == CigarOperator.X) return false;
		}
		return true;
	}
	
	private int getAnchorLength() {
		CigarElement ce;
		if (getBreakendDirection() == BreakendDirection.Forward) {
			ce = record.getCigar().getCigarElement(0);
		} else {
			ce = record.getCigar().getCigarElement(record.getCigarLength() - 1);
		}
		if (ce.getOperator() == CigarOperator.X) {
			return 0;
		}
		return record.getReadLength() - getBreakendLength();
	}
	private int getBreakendLength() {
		return getBreakendLength(record);
	}
	private static int getBreakendLength(SAMRecord record) {
		CigarElement ce;
		if (getBreakendDirection(record) == BreakendDirection.Forward) {
			ce = record.getCigar().getCigarElement(record.getCigarLength() - 1);
		} else {
			ce = record.getCigar().getCigarElement(0);
		}
		if (ce.getOperator() == CigarOperator.SOFT_CLIP) return ce.getLength();
		return 0;
	}
	private byte[] getAnchorBytes(byte[] data) {
		int len = getAnchorLength();
		if (getBreakendDirection() == BreakendDirection.Forward) {
			return Arrays.copyOfRange(data, 0, len);
		} else {
			return Arrays.copyOfRange(data, data.length - len, data.length);
		}
	}
	private byte[] getBreakendBytes(byte[] data) {
		int len = getBreakendLength();
		if (getBreakendDirection() == BreakendDirection.Forward) {
			return Arrays.copyOfRange(data, data.length - len, data.length);
		} else {
			return Arrays.copyOfRange(data, 0, len);
		}
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return breakend;
	}
	@Override
	public byte[] getBreakendSequence() {
		return getBreakendBytes(record.getReadBases());
	}
	@Override
	public byte[] getBreakendQuality() {
		return getBreakendBytes(record.getBaseQualities());
	}
	@Override
	public byte[] getAnchorSequence() {
		return getAssemblyAnchorSequence();
	}
	@Override
	public byte[] getAnchorQuality() {
		return getAssemblyAnchorQuals();
	}
	@Override
	public String getEvidenceID() {
		return record.getReadName();
	}
	@Override
	public AssemblyEvidenceSource getEvidenceSource() {
		return source;
	}
	@Override
	public int getLocalMapq() {
		return record.getMappingQuality();
	}
	@Override
	public byte[] getAssemblySequence() {
		if (isBreakendExact()) return record.getReadBases();
		// Need to remove placeholder Ns from inexact breakend
		return getBreakendDirection() == BreakendDirection.Forward ? Bytes.concat(getAssemblyAnchorSequence(), getBreakendSequence()) : Bytes.concat(getBreakendSequence(), getAssemblyAnchorSequence());
	}
	@Override
	public byte[] getAssemblyAnchorSequence() {
		return getAnchorBytes(record.getReadBases());
	}
	public byte[] getAssemblyAnchorQuals() {
		return getAnchorBytes(record.getBaseQualities());
	}
	@Override
	public int getAssemblyAnchorLength() {
		return getAnchorLength();
	}
	
	public int getAssemblySupportCount() {
		if (!isAnnotated()) {
			return getEvidenceIDs().size();
		} else {
			return getAssemblySupportCountSoftClip() + getAssemblySupportCountReadPair();
		}
	}
	@Override
	public int getAssemblySupportCountReadPair(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_READPAIR_COUNT), category, 0);
	}
	@Override
	public int getAssemblyReadPairLengthMax(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX), category, 0);
	}
	@Override
	public int getAssemblySupportCountSoftClip(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT), category, 0);
	}
	@Override
	public int getAssemblySupportCountSoftClipRemote(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_REMOTE_COUNT), category, 0);
	}
	public int getAssemblyNonSupportingReadPairCount(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT), category, 0);
	}
	public int getAssemblyNonSupportingReadPairCount() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblyNonSupportingSoftClipCount(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT), category, 0);
	}
	public int getAssemblyNonSupportingSoftClipCount() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT)).stream().mapToInt(x -> x).sum();
	}
	@Override
	public int getAssemblySoftClipLengthTotal(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL), category, 0);
	}
	@Override
	public int getAssemblySoftClipLengthMax(int category) {
		assertAnnotated();
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX), category, 0);
	}
	public float getAssemblySupportReadPairQualityScore(int category) {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_READPAIR_QUAL), category, 0);
	}
	public float getAssemblySupportSoftClipQualityScore(int category) {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL), category, 0);
	}
	public float getAssemblySupportRemoteQualityScore(int category) {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_REMOTE_QUAL), category, 0);
	}
	public float getAssemblyNonSupportingReadPairQualityScore(int category) {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL), category, 0);
	}
	public float getAssemblyNonSupportingReadPairQualityScore() {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblyNonSupportingSoftClipQualityScore(int category) {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL), category, 0);
	}
	public float getAssemblyNonSupportingSoftClipQualityScore() {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	@Override
	public int getAssemblySupportCountReadPair() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_READPAIR_COUNT)).stream().mapToInt(x -> x).sum();
	}
	@Override
	public int getAssemblyReadPairLengthMax() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX)).stream().mapToInt(x -> x).sum();
	}
	@Override
	public int getAssemblySupportCountSoftClip() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT)).stream().mapToInt(x -> x).sum();
	}
	@Override
	public int getAssemblySupportCountSoftClipRemote() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_REMOTE_COUNT)).stream().mapToInt(x -> x).sum();
	}
	@Override
	public int getAssemblySoftClipLengthTotal() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL)).stream().mapToInt(x -> x).sum();
	}
	@Override
	public int getAssemblySoftClipLengthMax() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX)).stream().mapToInt(x -> x).sum();
	}
	public float getAssemblySupportReadPairQualityScore() {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_READPAIR_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblySupportSoftClipQualityScore() {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblySupportSoftClipRemoteQualityScore() {
		assertAnnotated();
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_REMOTE_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblyNonSupportingQualityScore() {
		assertAnnotated();
		return (float)(AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL)).stream().mapToDouble(x -> x).sum() +
				AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL)).stream().mapToDouble(x -> x).sum());
	}
	public int getAssemblyNonSupportingCount() {
		assertAnnotated();
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT)).stream().mapToInt(x -> x).sum() +
				AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT)).stream().mapToInt(x -> x).sum();
	}
	private boolean isAnnotated() {
		Object attr = record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT);
		return attr != null;
	}
	private void assertAnnotated() {
		if (!isAnnotated()) {
			throw new IllegalStateException("Assembly has not yet been annotated");
		}
	}
	@Override
	public boolean isAssemblyFiltered() {
		return StringUtils.isNotBlank((String)record.getAttribute(SamTags.ASSEMBLY_FILTERS));
	}
	@Override
	public void filterAssembly(VcfFilter reason) {
		String tag = (String)record.getAttribute(SamTags.ASSEMBLY_FILTERS);
		if (StringUtils.isBlank(tag)) {
			tag = reason.filter();
		} else if (!tag.contains(reason.filter())) {
			tag = tag + "," + reason.filter();
		}
		getBackingRecord().setAttribute(SamTags.ASSEMBLY_FILTERS, tag);
		getSAMRecord().setAttribute(SamTags.ASSEMBLY_FILTERS, tag);
		getRemoteSAMRecord().setAttribute(SamTags.ASSEMBLY_FILTERS, tag);
	}
	@Override
	public List<VcfFilter> getFilters() {
		List<VcfFilter> list = Lists.newArrayList();
		String filters = (String)record.getAttribute(SamTags.ASSEMBLY_FILTERS);
		if (!StringUtils.isEmpty(filters)) {
			for (String s : filters.split(",")) {
				list.add(VcfFilter.get(s));
			}
		}
		return list;
	}
	public SAMRecord getSAMRecord() {
		return record;
	}
	public SAMRecord getRemoteSAMRecord() {
		return realignment;
	}
	public SAMRecord getBackingRecord() {
		return getSAMRecord();
	}
	@Override
	public String toString() {
		return String.format("A  %s N=%s", getBreakendSummary(), getEvidenceID());
	}
	@Override
	public boolean isBreakendExact() {
		return isExact;
	}
	/**
	 * Indicates that the breakend sequence was not aligned to the reference
	 * @return true if the breakend is likely to be novel sequence, false otherwise
	 */
	public boolean isNovelBreakend() {
		return realignment.getReadUnmappedFlag();
	}
	@Override
	public float getBreakendQual() {
		if (getBreakendLength() == 0) return 0;
		int rp = getAssemblySupportCountReadPair();
		double rpq = getAssemblySupportReadPairQualityScore();
		int sc = getAssemblySupportCountSoftClip();
		double scq =  getAssemblySupportSoftClipQualityScore();
		if (source.getContext().getAssemblyParameters().excludeNonSupportingEvidence) {
			rp -= getAssemblyNonSupportingReadPairCount();
			rpq -= getAssemblyNonSupportingReadPairQualityScore();
			sc -= getAssemblyNonSupportingSoftClipCount();
			scq -= getAssemblyNonSupportingSoftClipQualityScore();
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreBreakendAssembly(
				rp, rpq,
				sc, scq,
				getLocalMapq());
	}
	/**
	 * Performs local Smith-Watermann alignment of assembled contig to remove artifacts caused by misalignment of soft clipped reads.
	 * 
	 * Realignment of the full sequence (eg 50M100S) of small events results in unintended realignment
	 * to the far side of the event (eg 50S100M). The number of soft clipped bases aligned is capped to
	 * prevent this.
	 * 
	 * @param realignmentWindowSize number of additional reference bases to include around alignment
	 * @param portionOfAnchorRetained minimum portion of anchor bases that remain anchored to the reference after realignment
	 * @return includeBreakendBases include non-reference breakend bases when calculating alignment window size
	 */
	public SAMRecordAssemblyEvidence realign(int realignmentWindowSize, float minPortionOfAnchorRetained) {
		if (!isBreakendExact()) {
			return this;
		}
		SAMRecord read = getBackingRecord();
		ReferenceLookup reference = getEvidenceSource().getContext().getReference();
		SAMRecord realigned = SAMRecordUtil.realign(reference, read, realignmentWindowSize, true);
		if (!isValidRealignment(read.getCigar(), realigned.getCigar(), minPortionOfAnchorRetained)) {
			// Losing too many anchored bases can occur due to alignment to other side of small events
			// (since the breakend is usually longer than the anchor)
			// Try to just refine the breakend position by reducing the window size
        	realigned = SAMRecordUtil.realign(reference, read, realignmentWindowSize, false);
        	if (!isValidRealignment(read.getCigar(), realigned.getCigar(), minPortionOfAnchorRetained)) {
        		// can't find a valid realignment
        		return this;
        	}
        }
		if (realigned == getBackingRecord()) {
			// alread have the best alignment
			return this;
		}
		realigned.setReadName(getBackingRecord().getReadName() + "_r");
        SAMRecordAssemblyEvidence realignedAssembly = AssemblyFactory.hydrate(getEvidenceSource(), realigned);
        return realignedAssembly;
	}
	private static boolean isValidRealignment(Cigar original, Cigar realigned, float minPortionOfAnchorRetained) {
		return CigarUtil.commonReferenceBases(realigned, original) >= CigarUtil.countMappedBases(original.getCigarElements()) * minPortionOfAnchorRetained;
	}
	/**
	 * Determines whether the assembly is of the reference allele
	 * @return
	 */
	public boolean isReferenceAssembly() {
		return getBreakendLength() == 0 && CigarUtil.countIndels(getSAMRecord().getCigar()) == 0;
	}
}
