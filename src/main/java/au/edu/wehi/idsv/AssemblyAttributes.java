package au.edu.wehi.idsv;

import au.edu.wehi.idsv.debruijn.ContigCategorySupportHelper;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.function.DoubleBinaryOperator;
import java.util.function.IntBinaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AssemblyAttributes {
	private static final Log log = Log.getInstance(AssemblyAttributes.class);
	private static final String ID_COMPONENT_SEPARATOR = " ";
	private final SAMRecord record;
	private HashSet<String> evidenceIDs = null;
	private class Support {
		private final List<int[]> count;
		private final List<float[]> qual;
		public Support(SAMRecord record, String countAttr, String qualAttr) {
			int expectedLength = record.getReadLength() + 1;
			if (isUnanchored(record)) {
				expectedLength = Math.max(SAMRecordUtil.getStartClipLength(record), SAMRecordUtil.getEndClipLength(record)) + 1;
			}
			this.count = ContigCategorySupportHelper.unpackInts(record.getSignedIntArrayAttribute(countAttr), expectedLength);
			this.qual = ContigCategorySupportHelper.unpackFloats(record.getFloatArrayAttribute(qualAttr), expectedLength);
		}
	}
	private Support srSupport = null;
	private Support rpSupport = null;
	private Support getSplitReadSupport() {
		if (srSupport == null) {
			srSupport = new Support(record, SamTags.ASSEMBLY_SOFTCLIP_COUNT, SamTags.ASSEMBLY_SOFTCLIP_QUAL);
		}
		return srSupport;
	}
	private Support getReadPairSupport() {
		if (rpSupport == null) {
			rpSupport = new Support(record, SamTags.ASSEMBLY_READPAIR_COUNT, SamTags.ASSEMBLY_READPAIR_QUAL);
		}
		return rpSupport;
	}
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
			String[] ids = encoded.split(ID_COMPONENT_SEPARATOR);
			evidenceIDs = new HashSet<String>(Arrays.asList(ids));
			evidenceIDs.remove("");
		}
		return evidenceIDs;
	}
	public List<String> getOriginatingFragmentID() {
		String encoded = record.getStringAttribute(SamTags.ASSEMBLY_SUPPORTING_FRAGMENTS);
		if (encoded == null) {
			return ImmutableList.of();
		}
		return toList(encoded, ID_COMPONENT_SEPARATOR);
	}
	private static final List<String> toList(String str, String separator) {
		return Arrays.stream(str.split(ID_COMPONENT_SEPARATOR))
				.filter(s -> StringUtils.isNotEmpty(s))
				.collect(Collectors.toList());
	}
	/**
	 * Breakdown of DNA fragment support by category
	 * @param category
	 * @return
	 */
	public List<String> getOriginatingFragmentID(int category) {
		String encoded = record.getStringAttribute(SamTags.ASSEMBLY_SUPPORTING_FRAGMENTS);
		if (encoded == null) {
			return ImmutableList.of();
		}
		String[] categoryString = encoded.split(ID_COMPONENT_SEPARATOR + ID_COMPONENT_SEPARATOR);
		if (categoryString.length <= category) {
			return ImmutableList.of();
		}
		return Arrays.stream(categoryString[category].split(ID_COMPONENT_SEPARATOR))
				.filter(s -> StringUtils.isNotEmpty(s))
				.collect(Collectors.toList());
	}
	public static void annotatePerBasePerCategorySupport(
			SAMRecord record,
			List<int[]> scCount,
			List<float[]> scQual,
			List<int[]> rpCount,
			List<float[]> rpQual) {
		record.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_COUNT, ContigCategorySupportHelper.packInts(scCount));
		record.setAttribute(SamTags.ASSEMBLY_READPAIR_COUNT, ContigCategorySupportHelper.packInts(rpCount));
		record.setAttribute(SamTags.ASSEMBLY_SOFTCLIP_QUAL, ContigCategorySupportHelper.packFloats(scQual));
		record.setAttribute(SamTags.ASSEMBLY_READPAIR_QUAL, ContigCategorySupportHelper.packFloats(rpQual));
	}
	public static void annotatePerCategorySupport(
			SAMRecord record,
			List<Integer> scCount,
			List<Float> scQual,
			List<Integer> rpCount,
			List<Float> rpQual) {
		final int expectedLength = isUnanchored(record) ? Math.max(SAMRecordUtil.getStartClipLength(record), SAMRecordUtil.getEndClipLength(record)) + 1 : record.getReadLength() + 1;
		annotatePerBasePerCategorySupport(
			record,
			scCount.stream().map(x -> { int[] arr = new int[expectedLength]; Arrays.fill(arr, x); return arr; }).collect(Collectors.toList()),
			scQual.stream().map(x -> { float[] arr = new float[expectedLength]; Arrays.fill(arr, x); return arr; }).collect(Collectors.toList()),
			rpCount.stream().map(x -> { int[] arr = new int[expectedLength]; Arrays.fill(arr, x); return arr; }).collect(Collectors.toList()),
			rpQual.stream().map(x -> { float[] arr = new float[expectedLength]; Arrays.fill(arr, x); return arr; }).collect(Collectors.toList()));
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
		int maxLocalMapq = 0;
		int rpMaxLen = 0;
		List<Integer> scCount = Lists.newArrayList(Collections.nCopies(context.getCategoryCount(), 0));
		List<Float> scQual = Lists.newArrayList(Collections.nCopies(context.getCategoryCount(), 0f));
		List<Integer> rpCount = Lists.newArrayList(Collections.nCopies(context.getCategoryCount(), 0));
		List<Float> rpQual = Lists.newArrayList(Collections.nCopies(context.getCategoryCount(), 0f));
		for (DirectedEvidence e : support) {
			int category = ((SAMEvidenceSource)e.getEvidenceSource()).getSourceCategory();
			float qual = e.getBreakendQual();
			assert(e != null);
			maxLocalMapq = Math.max(maxLocalMapq, e.getLocalMapq());
			if (e instanceof NonReferenceReadPair) {
				rpMaxLen = Math.max(rpMaxLen, ((NonReferenceReadPair)e).getNonReferenceRead().getReadLength());
				rpCount.set(category, rpCount.get(category) + 1);
				rpQual.set(category, rpQual.get(category) + qual);
			} else {
				scCount.set(category, scCount.get(category) + 1);
				scQual.set(category, scQual.get(category) + qual);
			}
		}
		ensureUniqueEvidenceID(record.getReadName(), support);
		
		Map<Integer, List<DirectedEvidence>> evidenceByCategory = support.stream()
				.collect(Collectors.groupingBy(e -> ((SAMEvidenceSource)e.getEvidenceSource()).getSourceCategory()));
		String evidenceString = IntStream.range(0, evidenceByCategory.keySet().stream().mapToInt(x -> x).max().orElse(0) + 1)
			.mapToObj(i -> evidenceByCategory.get(i) == null ? "" : evidenceByCategory.get(i).stream()
					.map(e -> e.getEvidenceID())
					.distinct()
					.sorted()
					.collect(Collectors.joining(ID_COMPONENT_SEPARATOR)))
			.collect(Collectors.joining(ID_COMPONENT_SEPARATOR + ID_COMPONENT_SEPARATOR));
		String fragmentString = IntStream.range(0, evidenceByCategory.keySet().stream().mapToInt(x -> x).max().orElse(0) + 1)
			.mapToObj(i -> evidenceByCategory.get(i) == null ? "" : evidenceByCategory.get(i).stream()
					.flatMap(x -> x.getOriginatingFragmentID(i).stream())
					.distinct()
					.sorted()
					.collect(Collectors.joining(ID_COMPONENT_SEPARATOR)))
			.collect(Collectors.joining(ID_COMPONENT_SEPARATOR + ID_COMPONENT_SEPARATOR));
		record.setAttribute(SamTags.EVIDENCEID, evidenceString);
		record.setAttribute(SamTags.ASSEMBLY_SUPPORTING_FRAGMENTS, fragmentString);
		record.setAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX, rpMaxLen);
		record.setAttribute(SamTags.ASSEMBLY_STRAND_BIAS, (float)(support.size() == 0 ? 0.0 : (support.stream().mapToDouble(de -> de.getStrandBias()).sum() / support.size())));
		// TODO: proper mapq model
		record.setMappingQuality(maxLocalMapq);
		if (record.getMappingQuality() < context.getConfig().minMapq) {
			if (!MessageThrottler.Current.shouldSupress(log, "below minimum mapq")) {
				log.warn(String.format("Sanity check failure: %s has mapq below minimum", record.getReadName()));
			}
		}
		// pad out a default since we don't know any positional information here
		annotatePerCategorySupport(record, scCount, scQual, rpCount, rpQual);
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
	public enum SupportType {
		SplitRead,
		ReadPair,
	}
	public int getMinQualPosition(Range<Integer> assemblyContigOffset) {
		return getSupportingQualScore(assemblyContigOffset, null, null, Math::min).getLeft();
	}
	public Pair<Integer, Integer> getSupportingReadCount(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<SupportType> support, IntBinaryOperator positionalChoiceOperator) {
		List<List<int[]>> counts = new ArrayList<>();
		if (support == null || support.contains(SupportType.ReadPair)) {
			counts.add(getReadPairSupport().count);
		}
		if (support == null || support.contains(SupportType.SplitRead)) {
			counts.add(getSplitReadSupport().count);
		}
		int bestPos = -1;
		int best = -1;
		if (counts.size() == 0) {
			throw new IllegalArgumentException("support must be specified");
		}
		for (int i = assemblyContigOffset.lowerEndpoint(); i <= assemblyContigOffset.upperEndpoint(); i++) {
			int current = 0;
			for (int j = 0; j < counts.get(0).size(); j++) {
				if (supportingCategories == null || supportingCategories.contains(j)) {
					for (List<int[]> list : counts) {
						int[] arr = list.get(j);
						if (arr != null && arr.length > i) {
							current += arr[i];
						}
					}
				}
			}
			if (bestPos < 0) {
				bestPos = i;
				best = current;
			} else {
				int newBest = positionalChoiceOperator.applyAsInt(best, current);
				if (newBest != best) {
					bestPos = i;
					best = newBest;
				}
			}
		}
		return Pair.of(bestPos, best);
	}
	public Pair<Integer, Float> getSupportingQualScore(Range<Integer> assemblyContigOffset, List<Boolean> supportingCategories, Set<SupportType> support, DoubleBinaryOperator positionalChoiceOperator) {
		List<List<float[]>> quals = new ArrayList<>();
		if (support == null || support.contains(SupportType.ReadPair)) {
			quals.add(getReadPairSupport().qual);
		}
		if (support == null || support.contains(SupportType.SplitRead)) {
			quals.add(getSplitReadSupport().qual);
		}
		int bestPos = -1;
		float best = -1;
		if (quals.size() == 0) {
			throw new IllegalArgumentException("support must be specified");
		}
		for (int i = assemblyContigOffset.lowerEndpoint(); i <= assemblyContigOffset.upperEndpoint(); i++) {
			float current = 0;
			for (int j = 0; j < quals.get(0).size(); j++) {
				if (supportingCategories == null || supportingCategories.contains(j)) {
					for (List<float[]> list : quals) {
						float[] arr = list.get(j);
						if (arr != null && arr.length > i) {
							current += arr[i];
						}
					}
				}
			}
			if (bestPos < 0) {
				bestPos = i;
				best = current;
			} else {
				float newBest = (float)positionalChoiceOperator.applyAsDouble(best, current);
				if (newBest != best) {
					bestPos = i;
				}
			}
		}
		return Pair.of(bestPos, best);
	}
	public int getAssemblyReadPairLengthMax() {
		return record.getIntegerAttribute(SamTags.ASSEMBLY_READPAIR_LENGTH_MAX);
	}
	public BreakendDirection getAssemblyDirection() {
		Character c = (Character)record.getAttribute(SamTags.ASSEMBLY_DIRECTION);
		if (c == null) return null;
		return BreakendDirection.fromChar(c);
	}
	public double getStrandBias() {
		return AttributeConverter.asDouble(record.getAttribute(SamTags.ASSEMBLY_STRAND_BIAS), 0);
	}
}
