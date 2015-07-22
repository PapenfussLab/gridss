package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.sam.SamTags;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Sets;

/**
 * Tracks evidenceIDs backed by a SAMRecord
 * @author Daniel Cameron
 *
 */
public class EvidenceIDCollection {
	/**
	 * Using space as a separator as it is reserved in FASTA so shouldn't be in read names
	 */
	public static final String COMPONENT_EVIDENCEID_SEPARATOR = " ";
	private static final List<String> BACKING_ATTRIBUTES = ImmutableList.of(
			SamTags.ASSEMBLY_EVIDENCEID_PREFIX_DISCORDANT_PAIR,
			SamTags.ASSEMBLY_EVIDENCEID_PREFIX_UNMAPPED_MATE,
			SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ,
			SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE,
			SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SOFT_CLIP);
	/**
	 * Cache of evidence
	 */
	private List<HashMap<String, HashSet<String>>> byCategorySamTag = new ArrayList<HashMap<String,HashSet<String>>>();
	private HashSet<String> uncategorised;
	private HashSet<String> all = new HashSet<String>();
	public EvidenceIDCollection(int categoryCount, SAMRecord r) {
		for (int i = 0; i < categoryCount; i++) {		
			HashMap<String, HashSet<String>> lookup = new HashMap<String, HashSet<String>>();
			for (String samTag : BACKING_ATTRIBUTES) {
				HashSet<String> idSet = decodeEvidenceIDs(r.getAttribute(samTag + Integer.toString(i))); 
				lookup.put(samTag, idSet);
				all.addAll(idSet);
			}
			byCategorySamTag.add(lookup);
		}
		uncategorised = decodeEvidenceIDs(r.getAttribute(SamTags.ASSEMBLY_EVIDENCEID_UNCATEGORISED));
	}
	public EvidenceIDCollection(int categoryCount) {
		for (int i = 0; i < categoryCount; i++) {		
			HashMap<String, HashSet<String>> lookup = new HashMap<String, HashSet<String>>();
			for (String samTag : BACKING_ATTRIBUTES) {
				HashSet<String> idSet = decodeEvidenceIDs(null); 
				lookup.put(samTag, idSet);
				all.addAll(idSet);
			}
			byCategorySamTag.add(lookup);
		}
		uncategorised = new HashSet<String>();
	}
	public void add(String evidenceID) {
		uncategorised.add(evidenceID);
		all.add(evidenceID);
	}
	/**
	 * Adds the given evidence to the collection in the given category
	 * @param category ASSEMBLY_EVIDENCEID category defined in SamTags
	 * @param evidence evidenceID
	 */
	private void add(int categoryId, String type, String evidenceID) {		
		all.add(evidenceID);
		uncategorised.remove(evidenceID);
		byCategorySamTag.get(categoryId).get(type).add(evidenceID);
	}
	public Set<String> get() {
		return all;
	}
	/**
	 * Writes the evidence set to the given SAMRecord
	 */
	public void write(SAMRecord record) {
		record.setAttribute(SamTags.ASSEMBLY_EVIDENCEID_UNCATEGORISED, encodeEvidenceIDs(uncategorised));
		for (int i = 0; i < byCategorySamTag.size(); i++) {
			for (String backingPrefix : byCategorySamTag.get(i).keySet()) {
				String encoded = encodeEvidenceIDs(byCategorySamTag.get(i).get(backingPrefix));
				record.setAttribute(backingPrefix + Integer.toString(i), encoded);
			}
		}
	}
	private static String encodeEvidenceIDs(Collection<String> evidence) {
		List<String> evidenceIDs;
		if (evidence == null) {
			evidenceIDs = new ArrayList<String>();
		} else {
			evidenceIDs = new ArrayList<String>(evidence);
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
	private static HashSet<String> decodeEvidenceIDs(Object samTagAttributeValue) {
		HashSet<String> evidenceIds;
		if (samTagAttributeValue instanceof String && StringUtils.isNotEmpty((String)samTagAttributeValue)) {
			String[] ids = ((String)samTagAttributeValue).split(COMPONENT_EVIDENCEID_SEPARATOR);
			evidenceIds = Sets.newHashSet(ids);
		} else {
			evidenceIds = Sets.newHashSet();
		}
		return evidenceIds;
	}
	/**
	 * Categories the given evidence
	 * @param evidence
	 */
	public void categorise(DirectedEvidence evidence) {
		if (evidence == null) return;
		if (evidence instanceof AssemblyEvidence) return;
		int category = ((SAMEvidenceSource)evidence.getEvidenceSource()).getSourceCategory();
		if (evidence instanceof DiscordantReadPair) {
			add(category, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_DISCORDANT_PAIR, evidence.getEvidenceID());
		} else if (evidence instanceof UnmappedMateReadPair) {
			add(category, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_UNMAPPED_MATE, evidence.getEvidenceID());
		} else if (evidence instanceof RealignedRemoteSoftClipEvidence) {
			add(category, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE, evidence.getEvidenceID());
		} else if (evidence instanceof RealignedRemoteSoftClipEvidence) {
			add(category, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE, evidence.getEvidenceID());
		} else if (evidence instanceof RealignedSoftClipEvidence) {
			add(category, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ, evidence.getEvidenceID());
		} else if (evidence instanceof SoftClipEvidence) {
			add(category, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SOFT_CLIP, evidence.getEvidenceID());
		} else {
			throw new IllegalArgumentException(String.format("Unknown evidence type %s", evidence.getClass()));
		}
	}
	public void removeAll(Collection<String> evidenceID) {
		all.removeAll(evidenceID);
		uncategorised.removeAll(evidenceID);
		for (HashMap<String, HashSet<String>> lookup : byCategorySamTag) {
			for (HashSet<String> v : lookup.values()) {
				v.removeAll(evidenceID);
			}
		}
	}
	public int getCount(int category, String type) {
		return byCategorySamTag.get(category).get(type).size();
	}
	public void add(EvidenceIDCollection collection) {
		all.addAll(collection.all);
		uncategorised.addAll(collection.all);
		for (int i = 0; i < byCategorySamTag.size(); i++) {
			for (String backingPrefix : byCategorySamTag.get(i).keySet()) {
				byCategorySamTag.get(i).get(backingPrefix).addAll(collection.byCategorySamTag.get(i).get(backingPrefix));
			}
		}
	}
}
