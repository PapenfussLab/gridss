package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import htsjdk.samtools.SAMRecord;
import it.unimi.dsi.fastutil.objects.Object2FloatOpenHashMap;

/**
 * Tracks evidenceID and qualities
 * @author Daniel Cameron
 *
 */
public class EvidenceIDCollection {
	/**
	 * Using space as a separator as it is reserved in FASTA so shouldn't be in read names
	 */
	public static final String COMPONENT_EVIDENCEID_SEPARATOR = " ";
	private static final Map<String, List<String>> CATEGORY_ATTRIBUTES = ImmutableList.of(
			SamTags.EVIDENCEID_PREFIX_DISCORDANT_PAIR,
			SamTags.EVIDENCEID_PREFIX_UNMAPPED_MATE,
			SamTags.EVIDENCEID_PREFIX_SPLIT_READ,
			SamTags.EVIDENCEID_PREFIX_SPLIT_READ_REMOTE,
			SamTags.EVIDENCEID_PREFIX_SOFT_CLIP).stream()
			.collect(Collectors.toMap(prefix -> prefix,
					prefix -> IntStream.range(0, 10).mapToObj(i -> prefix + Integer.toString(i)).collect(Collectors.toList())));
	private static final List<String> BACKING_ATTRIBUTES = Stream.concat(ImmutableList.of(
			SamTags.EVIDENCEID_UNCATEGORISED,
			SamTags.EVIDENCEID_ASSEMBLY,
			SamTags.EVIDENCEID_ASSEMBLY_REMOTE,
			SamTags.EVIDENCEID_BREAKEND_ASSEMBLY)
			.stream(), CATEGORY_ATTRIBUTES.entrySet().stream().flatMap(e -> e.getValue().stream()))
			.collect(Collectors.toList());
	private static String getTag(int category, String samTag) {
		return CATEGORY_ATTRIBUTES.get(samTag).get(category);
	}
	/**
	 * Cache of evidence
	 */
	private HashMap<String, Object2FloatOpenHashMap<String>> bySamTag = new HashMap<String, Object2FloatOpenHashMap<String>>();
	private HashSet<String> all = new HashSet<String>();
	public EvidenceIDCollection(SAMRecord r) {
		for (String samTag : BACKING_ATTRIBUTES) {
			Object2FloatOpenHashMap<String> idSet = decodeEvidenceIDs(r.getAttribute(samTag));
			bySamTag.put(samTag, idSet);
			all.addAll(idSet.keySet());
		}
	}
	public EvidenceIDCollection() {
		for (String samTag : BACKING_ATTRIBUTES) {
			bySamTag.put(samTag, decodeEvidenceIDs(null));
		}
	}
	public void addUncategorised(String evidenceID) {
		add(SamTags.EVIDENCEID_UNCATEGORISED, evidenceID, 0);
	}
	/**
	 * Adds the given evidence to the collection in the given category
	 * @param f 
	 * @param category ASSEMBLY_EVIDENCEID category defined in SamTags
	 * @param evidence evidenceID
	 */
	private void add(String samTag, String evidenceID, float qual) {		
		all.add(evidenceID);
		if (!samTag.equals(SamTags.EVIDENCEID_UNCATEGORISED)) {
			bySamTag.get(SamTags.EVIDENCEID_UNCATEGORISED).remove(evidenceID);
		}
		bySamTag.get(samTag).put(evidenceID, qual);
	}
	private void add(String samTag, Object2FloatOpenHashMap<String> evidenceID, Function<String, String> idTransform) {	
		for (String id : evidenceID.keySet()) {
			add(samTag, idTransform.apply(id), evidenceID.getFloat(id));
		}
	}
	private void add(String samTag, Object2FloatOpenHashMap<String> evidenceID) {		
		all.addAll(evidenceID.keySet());
		if (!samTag.equals(SamTags.EVIDENCEID_UNCATEGORISED)) {
			bySamTag.get(SamTags.EVIDENCEID_UNCATEGORISED).keySet().removeAll(evidenceID.keySet());
		}
		bySamTag.get(samTag).putAll(evidenceID);
	}
	public Set<String> get() {
		return all;
	}
	/**
	 * Writes the evidence set to the given SAMRecord
	 */
	public void write(SAMRecord record) {
		for (String samTag : BACKING_ATTRIBUTES) {
			record.setAttribute(samTag, encodeEvidenceIDs(bySamTag.get(samTag)));
		}
	}
	private static String encodeEvidenceIDs(Object2FloatOpenHashMap<String> object2FloatOpenHashMap) {
		List<String> evidenceIDs;
		if (object2FloatOpenHashMap == null) {
			evidenceIDs = new ArrayList<String>();
		} else {
			evidenceIDs = new ArrayList<String>(object2FloatOpenHashMap.keySet());
		}
		// Set deterministic ordering of evidence to allow for SAM diff
		Collections.sort(evidenceIDs);
		StringBuilder sb = new StringBuilder();
		for (String id : evidenceIDs) {
			assert(!id.contains(COMPONENT_EVIDENCEID_SEPARATOR));
			sb.append(id);
			sb.append(COMPONENT_EVIDENCEID_SEPARATOR);
		}
		if (evidenceIDs.size() > 0) {
			sb.deleteCharAt(sb.length() - 1);
		}
		if (sb.toString().length() == 0) return null;
		return sb.toString();
	}
	private static Object2FloatOpenHashMap<String> decodeEvidenceIDs(Object samTagAttributeValue) {
		Object2FloatOpenHashMap<String> evidenceIds = new Object2FloatOpenHashMap<String>();
		if (samTagAttributeValue instanceof String && StringUtils.isNotEmpty((String)samTagAttributeValue)) {
			String[] ids = ((String)samTagAttributeValue).split(COMPONENT_EVIDENCEID_SEPARATOR);
			evidenceIds = new Object2FloatOpenHashMap<String>(ids.length);
			for (String s : ids) {
				evidenceIds.put(s, 0);
			}
		} else {
			evidenceIds = new Object2FloatOpenHashMap<String>();
		}
		return evidenceIds;
	}
	/**
	 * Categories the given evidence
	 * @param evidence
	 */
	public void categorise(DirectedEvidence evidence) {
		if (evidence == null) return;
		String id =  evidence.getEvidenceID();
		float qual = evidence.getBreakendQual();
		if (evidence instanceof DirectedBreakpoint) {
			qual = ((DirectedBreakpoint) evidence).getBreakpointQual();
		}
		String samTag = null;
		if (evidence instanceof AssemblyEvidence) {
			if (evidence instanceof DirectedBreakpoint) {
				if (evidence instanceof RemoteEvidence) {
					samTag = SamTags.EVIDENCEID_ASSEMBLY_REMOTE;
				} else {
					samTag = SamTags.EVIDENCEID_ASSEMBLY;
				}
			} else {
				samTag = SamTags.EVIDENCEID_BREAKEND_ASSEMBLY;
			}
		} else {
			int category = ((SAMEvidenceSource)evidence.getEvidenceSource()).getSourceCategory();
			if (evidence instanceof DiscordantReadPair) {
				samTag = getTag(category, SamTags.EVIDENCEID_PREFIX_DISCORDANT_PAIR);
			} else if (evidence instanceof UnmappedMateReadPair) {
				samTag = getTag(category, SamTags.EVIDENCEID_PREFIX_UNMAPPED_MATE);
			} else if (evidence instanceof RealignedRemoteSoftClipEvidence) {
				samTag = getTag(category, SamTags.EVIDENCEID_PREFIX_SPLIT_READ_REMOTE);
			} else if (evidence instanceof RealignedSoftClipEvidence) {
				samTag = getTag(category, SamTags.EVIDENCEID_PREFIX_SPLIT_READ);
			} else if (evidence instanceof SoftClipEvidence) {
				samTag = getTag(category, SamTags.EVIDENCEID_PREFIX_SOFT_CLIP);
			}
		}
		add(samTag, id, qual);
	}
	public void remove(String evidenceID) {
		all.remove(evidenceID);
		for (Object2FloatOpenHashMap<String> lookup : bySamTag.values()) {
			lookup.remove(evidenceID);
		}
	}
	public void add(EvidenceIDCollection collection) {
		all.addAll(collection.all);
		for (String samTag: bySamTag.keySet()) {
			bySamTag.get(samTag).putAll(collection.bySamTag.get(samTag));
		}
	}
	public void addRemote(EvidenceIDCollection e) {
		add(SamTags.EVIDENCEID_UNCATEGORISED, e.bySamTag.get(SamTags.EVIDENCEID_UNCATEGORISED));
		// breakend ids don't change
		add(SamTags.EVIDENCEID_BREAKEND_ASSEMBLY, e.bySamTag.get(SamTags.EVIDENCEID_BREAKEND_ASSEMBLY));
		for (String samTags : Iterables.concat(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SOFT_CLIP), CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_UNMAPPED_MATE))) {
			add(samTags, e.bySamTag.get(samTags));
		}
		// breakpoint ids swap with equivalent remote ids
		add(SamTags.EVIDENCEID_ASSEMBLY_REMOTE, e.bySamTag.get(SamTags.EVIDENCEID_ASSEMBLY), id -> "R" + id);
		add(SamTags.EVIDENCEID_ASSEMBLY, e.bySamTag.get(SamTags.EVIDENCEID_ASSEMBLY_REMOTE), id -> id.substring(1));
		for (String samTags : CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_DISCORDANT_PAIR)) {
			add(samTags, e.bySamTag.get(samTags), id -> swapPairSuffix(id));
		}
		for (int i = 0; i < 10; i++) {
			add(getTag(i, SamTags.EVIDENCEID_PREFIX_SPLIT_READ_REMOTE), e.bySamTag.get(getTag(i, SamTags.EVIDENCEID_PREFIX_SPLIT_READ)), id -> "R" + id);
			add(getTag(i, SamTags.EVIDENCEID_PREFIX_SPLIT_READ), e.bySamTag.get(getTag(i, SamTags.EVIDENCEID_PREFIX_SPLIT_READ_REMOTE)), id -> id.substring(1));
		}
	}
	private static String swapPairSuffix(String evidenceID) {
		if (evidenceID.endsWith(SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX)) {
			return evidenceID.substring(0, evidenceID.length() - SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX.length()) + SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX;
		} else if (evidenceID.endsWith(SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX)) {
			return evidenceID.substring(0, evidenceID.length() - SAMRecordUtil.SECOND_OF_PAIR_NAME_SUFFIX.length()) + SAMRecordUtil.FIRST_OF_PAIR_NAME_SUFFIX;
		}
		return evidenceID;
	}
	public float[] getQual(VcfAttributes attribute, int categoryCount) {
		switch (attribute) {
			case BREAKPOINT_ASSEMBLY_QUAL:
				return new float[] { sum(SamTags.EVIDENCEID_ASSEMBLY) };
			case BREAKPOINT_SPLITREAD_QUAL:
				return sum(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SPLIT_READ), categoryCount);
			case BREAKPOINT_READPAIR_QUAL:
				return sum(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_DISCORDANT_PAIR), categoryCount);
			case BREAKPOINT_ASSEMBLY_QUAL_REMOTE:
				return new float[] { sum(SamTags.EVIDENCEID_ASSEMBLY_REMOTE) };
			case BREAKPOINT_SPLITREAD_QUAL_REMOTE:
				return sum(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SPLIT_READ_REMOTE), categoryCount);
			case BREAKEND_ASSEMBLY_QUAL:
				return new float[] { sum(SamTags.EVIDENCEID_BREAKEND_ASSEMBLY) };
			case BREAKEND_SOFTCLIP_QUAL:
				return sum(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SOFT_CLIP), categoryCount);
			case BREAKEND_UNMAPPEDMATE_QUAL:
				return sum(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_UNMAPPED_MATE), categoryCount);
		default:
			throw new IllegalArgumentException();
		}
	}
	private float[] sum(List<String> samTag, int categoryCount) {
		float[] result = new float[categoryCount];
		for (int i = 0; i <categoryCount; i++) {
			result[i] = sum(samTag.get(i));
		}
		return result;
	}
	private float sum(String samTag) {
		float sum = 0;
		for (float f : bySamTag.get(samTag).values()) sum += f;
		return sum;
	}
	private int[] count(List<String> samTag, int categoryCount) {
		int[] result = new int[categoryCount];
		for (int i = 0; i < categoryCount; i++) {
			result[i] = count(samTag.get(i));
		}
		return result;
	}
	private int count(String samTag) {
		return bySamTag.get(samTag).size();
	}
	public int[] getCount(VcfAttributes attribute, int categoryCount) {
		switch (attribute) {
			case BREAKPOINT_ASSEMBLY_COUNT:
				return new int[] { count(SamTags.EVIDENCEID_ASSEMBLY) };
			case BREAKPOINT_SPLITREAD_COUNT:
				return count(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SPLIT_READ), categoryCount);
			case BREAKPOINT_READPAIR_COUNT:
				return count(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_DISCORDANT_PAIR), categoryCount);
			case BREAKPOINT_ASSEMBLY_COUNT_REMOTE:
				return new int[] { count(SamTags.EVIDENCEID_ASSEMBLY_REMOTE) };
			case BREAKPOINT_SPLITREAD_COUNT_REMOTE:
				return count(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SPLIT_READ_REMOTE), categoryCount);
			case BREAKEND_ASSEMBLY_COUNT:
				return new int[] { count(SamTags.EVIDENCEID_BREAKEND_ASSEMBLY) };
			case BREAKEND_SOFTCLIP_COUNT:
				return count(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_SOFT_CLIP), categoryCount);
			case BREAKEND_UNMAPPEDMATE_COUNT:
				return count(CATEGORY_ATTRIBUTES.get(SamTags.EVIDENCEID_PREFIX_UNMAPPED_MATE), categoryCount);
		default:
			throw new IllegalArgumentException();
		}
	}
}
