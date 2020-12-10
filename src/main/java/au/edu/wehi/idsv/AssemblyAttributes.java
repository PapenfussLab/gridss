package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.Range;
import com.google.common.collect.Streams;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AssemblyAttributes {
	private static final Log log = Log.getInstance(AssemblyAttributes.class);
	private static final String ID_COMPONENT_SEPARATOR = " ";
	private final SAMRecord record;
	private Collection<AssemblyEvidenceSupport> support = null;
	public static boolean isAssembly(SAMRecord record) {
		return record.hasAttribute(SamTags.IS_ASSEMBLY);
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
		if (!isAssembly(record)) {
			throw new IllegalArgumentException("record is not an assembly.");
		}
		this.record = record;
	}
	public AssemblyAttributes(SingleReadEvidence record) {
		this(record.getSAMRecord());
	}

	public static void adjustAssemblyAnnotationDueToContigChange(SAMRecord record, int startTruncatedBases) {
		int[] intervalStart = record.getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START);
		int[] intervalEnd = record.getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END);
		if (intervalStart != null && intervalEnd != null) {
			for (int i = 0; i < intervalStart.length; i++) {
				intervalStart[i] -= startTruncatedBases;
				intervalEnd[i] -= startTruncatedBases;
			}
			record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, intervalStart);
			record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, intervalEnd);
		}
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
		return getSupport(e.getEvidenceSource() instanceof AssemblyEvidenceSource ? (AssemblyEvidenceSource)e.getEvidenceSource() : null)
				.stream()
				.anyMatch(ee -> ee.getEvidenceID().equals(e.getEvidenceID()));
	}
	private Collection<AssemblyEvidenceSupport> getSupport(AssemblyEvidenceSource aes) {
		if (support == null) {
			support = new ArrayList<>();
			if (!record.hasAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY)) {
				// can't write zero length SAM arrays
				// no attribute means no supporting evidence
			} else {
				byte[] type = record.getSignedByteArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE);
				int[] category = record.getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY);
				int[] intervalStart = record.getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START);
				int[] intervalEnd = record.getSignedIntArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END);
				float[] qual = record.getFloatArrayAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL);
				String[] evidenceId = record.getStringAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID).split(ID_COMPONENT_SEPARATOR);
				String[] fragmentId = record.getStringAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID).split(ID_COMPONENT_SEPARATOR);
				if (type == null || category == null || intervalStart == null || intervalEnd == null || qual == null) {
					String msg = "Sanity check failure:" + record.getReadName() + " missing required evidence SAM tag.";
					log.error(msg);
					throw new IllegalStateException(msg);
				}
				if (aes != null) {
					final int[] lookup = aes.getAssemblyCategoryToProcessingContextCategoryLookup();
					for (int i = 0; i < category.length; i++) {
						category[i] = lookup[category[i]];
					}
				}
				if (type.length != category.length
						|| type.length != intervalStart.length
						|| type.length != intervalEnd.length
						|| type.length != qual.length
						|| type.length != evidenceId.length
						|| type.length != fragmentId.length) {
					String msg = "Sanity check failure:" + record.getReadName() + " has inconsistent evidence SAM tag.";
					log.error(msg);
					throw new IllegalStateException(msg);
				}
				for (int i = 0; i < category.length; i++) {
					support.add(new AssemblyEvidenceSupport(
							AssemblyEvidenceSupport.SupportType.value(type[i]),
							Range.closed(intervalStart[i], intervalEnd[i]),
							evidenceId[i],
							fragmentId[i],
							category[i],
							qual[i]
					));
				}
			}
		}
		return support;
	}
	private static int maxReadLength(Collection<DirectedEvidence> support) {
		return support.stream()
				.mapToInt(e ->  e instanceof NonReferenceReadPair ? ((NonReferenceReadPair)e).getNonReferenceRead().getReadLength() : ((SingleReadEvidence)e).getSAMRecord().getReadLength())
				.max()
				.orElse(0);
	}
	private static float readStrandBias(Collection<DirectedEvidence> fullSupport) {
		List<DirectedEvidence> support = fullSupport.stream()
				.filter(s -> s instanceof SingleReadEvidence)
				.collect(Collectors.toList());
		if (support.size() == 0) {
			return 0.5f;
		}
		return (float)support.stream()
				.mapToDouble(e -> e.getStrandBias())
				.sum() / support.size();
	}
	private static int maxLocalMapq(Collection<DirectedEvidence> support) {
		return support.stream()
				.mapToInt(e -> e.getLocalMapq())
				.max()
				.orElse(0);
	}
	/**
	 * Annotates an assembly with summary information regarding the reads used to produce the assembly
	 */
	public static void annotateAssembly(ProcessingContext context, SAMRecord record, List<DirectedEvidence> support, List<AssemblyEvidenceSupport> aes) {
		if (support == null) {
			if (!MessageThrottler.Current.shouldSupress(log, "assemblies with no support")) {
				log.error("No support for assembly " + record.getReadName());
			}
			support = Collections.emptyList();
		}
		// #260 support and aes do not have to match as aes is filtered based on supporting interval
		// A read supporting a kmer but not the contig will have no aes support record
		//if (support.size() != aes.size()) {
		//	throw new IllegalArgumentException("support and aes sizes do not match");
		//}
		annotateAssemblyEvidenceSupport(record, aes);

		record.setAttribute(SamTags.IS_ASSEMBLY, (byte)1);
		record.setAttribute(SamTags.ASSEMBLY_MAX_READ_LENGTH, maxReadLength(support));
		record.setAttribute(SamTags.ASSEMBLY_STRAND_BIAS, readStrandBias(support));
		ensureUniqueEvidenceID(record.getReadName(), support);
		// TODO: proper mapq model
		record.setMappingQuality(maxLocalMapq(support));
		if (record.getMappingQuality() < context.getConfig().minMapq) {
			if (!MessageThrottler.Current.shouldSupress(log, "below minimum mapq")) {
				log.warn(String.format("Sanity check failure: %s has mapq below minimum", record.getReadName()));
			}
		}
	}

	private static void annotateAssemblyEvidenceSupport(SAMRecord record, List<AssemblyEvidenceSupport> aes) {
		if (aes == null || aes.size() == 0) {
			return;
		}
		Collections.sort(aes, AssemblyEvidenceSupport.ByEvidenceID);
		byte[] type = new byte[aes.size()];
		int[] category = new int[aes.size()];
		int[] intervalStart = new int[aes.size()];
		int[] intervalEnd = new int[aes.size()];
		float[] qual = new float[aes.size()];
		StringBuilder evidenceId = new StringBuilder();
		StringBuilder fragmentId = new StringBuilder();
		for (int i = 0; i < aes.size(); i++) {
			AssemblyEvidenceSupport s = aes.get(i);
			type[i] = (byte)s.getSupportType().getValue();
			category[i] = s.getCategory();
			intervalStart[i] = s.getAssemblyContigOffset().lowerEndpoint();
			intervalEnd[i] = s.getAssemblyContigOffset().upperEndpoint();
			qual[i] = s.getQual();
			if (i != 0) {
				evidenceId.append(ID_COMPONENT_SEPARATOR);
				fragmentId.append(ID_COMPONENT_SEPARATOR);
			}
			evidenceId.append(s.getEvidenceID());
			fragmentId.append(s.getFragmentID());
		}
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, type);
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, category);
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, evidenceId.toString());
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, fragmentId.toString());
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, intervalStart);
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, intervalEnd);
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, qual);
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
	private Stream<AssemblyEvidenceSupport> filterSupport(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		Stream<AssemblyEvidenceSupport> stream = getSupport(aes).stream();
		if (assemblyContigOffset != null) {
			stream = stream.filter(s -> s.getAssemblyContigOffset().isConnected(assemblyContigOffset));
		}
		if (supportingCategories != null) {
			stream = stream.filter(s -> supportingCategories.contains(s.getCategory()));
		}
		if (supportTypes != null) {
			stream = stream.filter(s -> supportTypes.contains(s.getSupportType()));
		}
		return stream;
	}
	public Collection<String> getEvidenceIDs(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		return filterSupport(assemblyContigOffset, supportingCategories, supportTypes, aes)
				.map(s -> s.getEvidenceID())
				.collect(Collectors.toList());
	}
	public Set<String> getOriginatingFragmentID(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		return filterSupport(assemblyContigOffset, supportingCategories, supportTypes, aes)
				.map(s -> s.getFragmentID())
				.collect(Collectors.toSet());
	}
	public int getMinQualPosition(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		if (assemblyContigOffset == null) {
			throw new NullPointerException("assemblyContigOffset is required.");
		}
		float best = getSupportingQualScore(assemblyContigOffset.lowerEndpoint(), supportingCategories, supportTypes, aes);
		int bestPos = assemblyContigOffset.lowerEndpoint();
		for (int i = assemblyContigOffset.lowerEndpoint() + 1; i <= assemblyContigOffset.upperEndpoint(); i++) {
			float current = getSupportingQualScore(i, null, null, aes);
			if (current < best) {
				best = current;
				bestPos = i;
			}
		}
		return bestPos;
	}
	public int getMaxQualPosition(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		if (assemblyContigOffset == null) {
			throw new NullPointerException("assemblyContigOffset is required.");
		}
		float best = getSupportingQualScore(assemblyContigOffset.lowerEndpoint(), supportingCategories, supportTypes, aes);
		int bestPos = assemblyContigOffset.lowerEndpoint();
		for (int i = assemblyContigOffset.lowerEndpoint() + 1; i <= assemblyContigOffset.upperEndpoint(); i++) {
			float current = getSupportingQualScore(i, null, null, aes);
			if (current > best) {
				best = current;
				bestPos = i;
			}
		}
		return bestPos;
	}
	public int getSupportingReadCount(Range<Integer> assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		return (int)filterSupport(assemblyContigOffset, supportingCategories, supportTypes, aes).count();
	}
	public int getSupportingReadCount(int assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		return (int)filterSupport(Range.closed(assemblyContigOffset, assemblyContigOffset), supportingCategories, supportTypes, aes).count();
	}
	public float getSupportingQualScore(int assemblyContigOffset, Set<Integer> supportingCategories, Set<AssemblyEvidenceSupport.SupportType> supportTypes, AssemblyEvidenceSource aes) {
		return (float)filterSupport(Range.closed(assemblyContigOffset, assemblyContigOffset), supportingCategories, supportTypes, aes).mapToDouble(s -> s.getQual()).sum();
	}
	public int getAssemblyMaxReadLength() {
		return record.getIntegerAttribute(SamTags.ASSEMBLY_MAX_READ_LENGTH);
	}
	public BreakendDirection getAssemblyDirection() {
		Character c = (Character)record.getAttribute(SamTags.ASSEMBLY_DIRECTION);
		if (c == null) return null;
		return BreakendDirection.fromChar(c);
	}
	public double getStrandBias() {
		return AttributeConverter.asDouble(record.getAttribute(SamTags.ASSEMBLY_STRAND_BIAS), 0);
	}
	public static int getUnanchoredPlacholderAnchoredBases(SAMRecord record) {
		if (!isUnanchored(record)) {
			throw new IllegalArgumentException("Assembly is not unanchored.");
		}
		Cigar c = getUnanchoredCigar(record);
		int xbases = c.getCigarElements().stream()
				.filter(ce -> ce.getOperator() == CigarOperator.X)
				.mapToInt(ce -> ce.getLength())
				.sum();
		return xbases;
	}
	private static Cigar getUnanchoredCigar(SAMRecord record) {
		return Streams.concat(Stream.of(new ChimericAlignment(record)),
				ChimericAlignment.getChimericAlignments(record).stream())
				.filter(ca -> ca.cigar.getCigarElements().stream().anyMatch(ce -> ce.getOperator() == CigarOperator.X))
				.map(ca -> ca.cigar)
				.findFirst()
				.orElse(null);
	}
}
