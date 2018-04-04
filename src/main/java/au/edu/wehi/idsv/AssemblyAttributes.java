package au.edu.wehi.idsv;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.NotImplementedException;

import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.MessageThrottler;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

public class AssemblyAttributes {
	private static final Log log = Log.getInstance(AssemblyAttributes.class);
	private static final String COMPONENT_EVIDENCEID_SEPARATOR = " ";
	private final SAMRecord record;
	private HashSet<String> evidenceIDs = null;
	public static boolean isAssembly(SAMRecord record) {
		return record.getAttribute(SamTags.EVIDENCEID) != null;
	}
	public static boolean isUnanchored(SAMRecord record) {
		return record.hasAttribute(SamTags.UNANCHORED);
		//return Iterables.any(record.getCigar().getCigarElements(), ce -> ce.getOperator() == CigarOperator.X);
	}
	public static boolean isAssembly(DirectedEvidence record) {
		if (record instanceof SingleReadEvidence) {
			return isAssembly(((SingleReadEvidence)record).getSAMRecord());
		}
		return false;
	}
	public AssemblyAttributes(SAMRecord record) {
		this.record = record;
	}
	public AssemblyAttributes(SingleReadEvidence record) {
		this(record.getSAMRecord());
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
	public Collection<String> getEvidenceIDs() {
		if (evidenceIDs == null) {
			String encoded = record.getStringAttribute(SamTags.EVIDENCEID);
			if (encoded == null) {
				throw new IllegalStateException("Unable to get constituent evidenceIDs from assembly with evidence tracking disabled");
			}
			String[] ids = encoded.split(COMPONENT_EVIDENCEID_SEPARATOR);
			evidenceIDs = new HashSet<String>(Arrays.asList(ids));
			evidenceIDs.remove("");
		}
		return evidenceIDs;
	}
	public static void annotateNonSupporting(ProcessingContext context, BreakpointSummary assemblyBreakpoint, SAMRecord record, Collection<DirectedEvidence> support) {
		int n = context.getCategoryCount();
		float[] nsrpQual = new float[n];
		float[] nsscQual = new float[n];
		int[] nsrpCount = new int[n];
		int[] nsscCount = new int[n];
		BreakpointSummary breakendWithMargin = (BreakpointSummary)context.getVariantCallingParameters().withMargin(assemblyBreakpoint);
		for (DirectedEvidence e : support) {
			int offset = ((SAMEvidenceSource)e.getEvidenceSource()).getSourceCategory();
			float qual = e.getBreakendQual();
			if (e instanceof NonReferenceReadPair) {
				if (breakendWithMargin != null && !breakendWithMargin.overlaps(e.getBreakendSummary())) {
					nsrpCount[offset]++;
					nsrpQual[offset] += qual;
				}
			} else if (e instanceof SingleReadEvidence) {
				if (breakendWithMargin != null && !breakendWithMargin.overlaps(e.getBreakendSummary())) {
					nsscCount[offset]++;
					nsscQual[offset] += qual;
				}
			} else {
				throw new NotImplementedException("Sanity check failure: not a read or a read pair.");
			}
		}
		record.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT, nsrpCount);
		record.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT, nsscCount);
		record.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL, nsrpQual);
		record.setAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL, nsscQual);
	}
	/**
	 * Annotates an assembly with summary information regarding the reads used to produce the assembly
	 */
	public static void annotateAssembly(ProcessingContext context, SAMRecord record, Collection<DirectedEvidence> support) {
		if (support == null) {
			if (!MessageThrottler.Current.shouldSupress(log, "assemblies with no support")) {
				log.error("No support for assembly " + record.getReadName());
			}
			support = Collections.emptyList();
		}
		int n = context.getCategoryCount();
		float[] rpQual = new float[n];
		float[] scQual = new float[n];
		int[] rpCount = new int[n];
		int[] rpMaxLen = new int[n];
		int[] scCount = new int[n];
		int[] scLenMax = new int[n];
		int[] scLenTotal = new int[n];
		int maxLocalMapq = 0;
		for (DirectedEvidence e : support) {
			assert(e != null);
			maxLocalMapq = Math.max(maxLocalMapq, e.getLocalMapq());
			int offset = ((SAMEvidenceSource)e.getEvidenceSource()).getSourceCategory();
			float qual = e.getBreakendQual();
			if (e instanceof NonReferenceReadPair) {
				rpCount[offset]++;
				rpQual[offset] += qual;
				rpMaxLen[offset] = Math.max(rpMaxLen[offset], ((NonReferenceReadPair)e).getNonReferenceRead().getReadLength());
			} else if (e instanceof SingleReadEvidence) {
				scCount[offset]++;
				scQual[offset] += qual;
				int clipLength = e.getBreakendSequence().length;
				scLenMax[offset] = Math.max(scLenMax[offset], clipLength);
				scLenTotal[offset] += clipLength;
			} else {
				throw new NotImplementedException("Sanity check failure: not a read or a read pair.");
			}
		}
		ensureUniqueEvidenceID(record.getReadName(), support);
		String evidenceString = support.stream()
				.map(e -> e.getEvidenceID())
				.sorted()
				.collect(Collectors.joining(COMPONENT_EVIDENCEID_SEPARATOR));
		record.setAttribute(SamTags.EVIDENCEID, evidenceString);
		record.setAttribute(SamTags.ASSEMBLY_READPAIR_COUNT, rpCount);
		record.setAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX, rpMaxLen);
		record.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT, scCount);
		record.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, scLenMax);
		record.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, scLenTotal);
		record.setAttribute(SamTags.ASSEMBLY_READPAIR_QUAL, rpQual);
		record.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL, scQual);
		// TODO: proper mapq model
		record.setMappingQuality(maxLocalMapq);
		if (record.getMappingQuality() < context.getConfig().minMapq) {
			if (!MessageThrottler.Current.shouldSupress(log, "below minimum mapq")) {
				log.warn(String.format("Sanity check failure: %s has mapq below minimum", record.getReadName()));
			}
		}
	}
	private static boolean ensureUniqueEvidenceID(String assemblyName, Collection<DirectedEvidence> support) {
		boolean isUnique = true;
		Set<String> map = new HashSet<String>();
		for (DirectedEvidence id : support) {
			if (map.contains(id.getEvidenceID())) {
				if (!MessageThrottler.Current.shouldSupress(log, "duplicated evidenceIDs")) {
					log.error("Found evidenceID " + id.getEvidenceID() + " multiple times in assembly " + assemblyName);
				}
				isUnique = false;
			}
			map.add(id.getEvidenceID());
		}
		return isUnique;
	}
	public int getAssemblySupportCount() {
		return getAssemblySupportCountSoftClip() + getAssemblySupportCountReadPair();
	}
	public int getAssemblySupportCountReadPair(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_READPAIR_COUNT), category, 0);
	}
	public int getAssemblyReadPairLengthMax(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX), category, 0);
	}
	public int getAssemblySupportCountSoftClip(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT), category, 0);
	}
	public int getAssemblyNonSupportingReadPairCount(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT), category, 0);
	}
	public int getAssemblyNonSupportingReadPairCount() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblyNonSupportingSoftClipCount(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT), category, 0);
	}
	public int getAssemblyNonSupportingSoftClipCount() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblySoftClipLengthTotal(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL), category, 0);
	}
	public int getAssemblySoftClipLengthMax(int category) {
		return AttributeConverter.asIntListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX), category, 0);
	}
	public float getAssemblySupportReadPairQualityScore(int category) {
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_READPAIR_QUAL), category, 0);
	}
	public float getAssemblySupportSoftClipQualityScore(int category) {
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL), category, 0);
	}
	public float getAssemblyNonSupportingReadPairQualityScore(int category) {
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL), category, 0);
	}
	public float getAssemblyNonSupportingReadPairQualityScore() {
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblyNonSupportingSoftClipQualityScore(int category) {
		return (float)AttributeConverter.asDoubleListOffset(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL), category, 0);
	}
	public float getAssemblyNonSupportingSoftClipQualityScore() {
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public int getAssemblySupportCountReadPair() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_READPAIR_COUNT)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblyReadPairLengthMax() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblySupportCountSoftClip() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblySoftClipLengthTotal() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL)).stream().mapToInt(x -> x).sum();
	}
	public int getAssemblySoftClipLengthMax() {
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX)).stream().mapToInt(x -> x).sum();
	}
	public float getAssemblySupportReadPairQualityScore() {
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_READPAIR_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblySupportSoftClipQualityScore() {
		return (float)AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL)).stream().mapToDouble(x -> x).sum();
	}
	public float getAssemblyNonSupportingQualityScore() {
		return (float)(AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_QUAL)).stream().mapToDouble(x -> x).sum() +
				AttributeConverter.asDoubleList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_QUAL)).stream().mapToDouble(x -> x).sum());
	}
	public int getAssemblyNonSupportingCount() {		
		return AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_READPAIR_COUNT)).stream().mapToInt(x -> x).sum() +
				AttributeConverter.asIntList(record.getAttribute(SamTags.ASSEMBLY_NONSUPPORTING_SOFTCLIP_COUNT)).stream().mapToInt(x -> x).sum();
	}
	public BreakendDirection getAssemblyDirection() {
		Character c = (Character)record.getAttribute(SamTags.ASSEMBLY_DIRECTION);
		if (c == null) return null;
		return BreakendDirection.fromChar((char)c);
	}
}
