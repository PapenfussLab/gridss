package au.edu.wehi.idsv;

import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;

public class GreedyVariantAllocationCache {
	/**
	 * read (pair) -> (event, score, read pair alignment)
	 * Only best placement of the read pairs should be allocated.
	 */
	private final HashMap<Hash128bit, Node> bestReadPairAlignment;
	/**
	 * read -> (event, score, read alignment)
	 * 
	 * Only evidence from the best alignment for each read should be allocated
	 * since we don't want to allocate the mutually exclusive evidence that
	 * results from two separate read alignments for the same read
	 */
	private final HashMap<Hash128bit, Node> bestReadAlignment;
	/**
	 * evidenceID -> best event lookup
	 * Each evidence can support only a single variant.
	 */
	private final HashMap<Hash128bit, Node> bestEventForEvidence;
	public GreedyVariantAllocationCache(boolean ensureUniqueReadPairAlignment, boolean ensureUniqueReadAlignment, boolean ensureUniqueEvidenceAllocation) {
		this.bestReadPairAlignment = ensureUniqueReadPairAlignment ? new HashMap<>() : null;
		this.bestReadAlignment = ensureUniqueReadAlignment ? new HashMap<>() : null;
		this.bestEventForEvidence = ensureUniqueEvidenceAllocation ? new HashMap<>() : null;
	}
	public void addBreakpoint(String event, float score, DirectedEvidence evidence) {
		addBreakpoint(new Hash128bit(event), score, evidence);
	}
	private void addBreakpoint(Hash128bit event, float score, DirectedEvidence evidence) {
		put(bestEventForEvidence, new Hash128bit(evidence.getEvidenceID()), null, event, score);
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			Hash128bit readpairid = new Hash128bit(dp.getLocalledMappedRead().getReadName());
			Hash128bit alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			put(bestReadPairAlignment, readpairid, alignment, event, score);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			Hash128bit readid = new Hash128bit(r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r)));
			Hash128bit alignment = new Hash128bit(getReadAlignment(r));
			put(bestReadAlignment, readid, alignment, event, score);
		}
	}
	private String getReadAlignment(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return "";
		List<ChimericAlignment> ca = Lists.newArrayList(ChimericAlignment.getChimericAlignments(r));
		ca.add(new ChimericAlignment(r));
		return getAlignmentString(ca);
	}
	private Hash128bit getReadPairAlignment(SAMRecord r) {
		List<ChimericAlignment> ca = new ArrayList<>(2);
		if (!r.getReadUnmappedFlag()) {
			ca.add(new ChimericAlignment(r));
		}
		if (!r.getMateUnmappedFlag()) {
			ca.add(new ChimericAlignment(
					r.getMateReferenceName(),
					r.getMateAlignmentStart(),
					r.getMateNegativeStrandFlag(),
					TextCigarCodec.decode(r.getStringAttribute(SAMTag.MC.name())),
					// not using these fields so zeros are fine
					0,
					0));
		}
		return new Hash128bit(getAlignmentString(ca));
	}
	private static final Comparator<ChimericAlignment> ByPosition = Comparator.<ChimericAlignment, Integer>comparing(ca -> ca.pos)
			.thenComparing(ca -> ca.isNegativeStrand)
			.thenComparing(ca -> ca.rname)
			.thenComparing(ca -> ca.cigar.toString());
	private String getAlignmentString(List<ChimericAlignment> ca) {
		StringBuilder sb = new StringBuilder();
		// need to make sure the order is the same at every alignment
		Collections.sort(ca, ByPosition);
		for (ChimericAlignment a : ca) {
			sb.append(a.rname);
			sb.append('#');
			sb.append(a.isNegativeStrand ? '-' : '+');
			sb.append(a.pos);
			sb.append('#');
			sb.append(a.cigar.toString());
			sb.append('#');
		}
		return sb.toString();
	}
	private static Hash128bit getEvent(VariantContextDirectedBreakpoint variant) {
		return new Hash128bit(variant.getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, null));
	}
	public void addBreakpoint(VariantContextDirectedBreakpoint variant, List<DirectedEvidence> evidence) {
		Hash128bit event = getEvent(variant);
		for (DirectedEvidence e : evidence) {
			addBreakpoint(event, variant.getBreakpointQual(), e);
		}
	}
	private void put(HashMap<Hash128bit, Node> bestAlignment, Hash128bit readName, Hash128bit alignment, Hash128bit event, float score) {
		if (bestAlignment == null) return;
		Node existing = bestAlignment.get(readName);
		if (existing == null || existing.score < score) {
			bestAlignment.put(readName, new Node(event, score, alignment));
		}
	}
	private boolean isBestAlignment(HashMap<Hash128bit, Node> bestAlignment, Hash128bit readName, Hash128bit alignment) {
		if (bestAlignment == null) return true;
		Node node = bestAlignment.get(readName);
		return node != null && Objects.equals(alignment, node.alignment);
	}
	public boolean isBestBreakpoint(VariantContextDirectedBreakpoint variant, DirectedEvidence evidence) {
		return isBestBreakpoint(getEvent(variant), evidence);
	}
	public boolean isBestBreakpoint(String event, DirectedEvidence evidence) {
		return isBestBreakpoint(new Hash128bit(event), evidence);
	}
	public boolean isBestBreakpoint(Hash128bit event, DirectedEvidence evidence) {
		if (bestEventForEvidence != null && !event.equals(bestEventForEvidence.get(new Hash128bit(evidence.getEvidenceID())).event)) {
			// This is not the best breakpoint supported by this evidence
			return false;
		}
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			Hash128bit readpairid = new Hash128bit(dp.getLocalledMappedRead().getReadName());
			Hash128bit alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			return isBestAlignment(bestReadPairAlignment, readpairid, alignment);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			Hash128bit readid = new Hash128bit(r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r)));
			Hash128bit alignment = new Hash128bit(getReadAlignment(r));
			return isBestAlignment(bestReadAlignment, readid, alignment);
		}
	}
	private static class Node {
		public Node(Hash128bit event, float score, Hash128bit alignment) {
			this.event = event;
			this.score = score;
			this.alignment = alignment;
		}
		public final Hash128bit event;
		public final Hash128bit alignment;
		public final float score;
		public String toString() {
			return String.format("(%s,%f,%s)", event, score, alignment);
		}
	}
	private static class Hash128bit {
		private final long key1;
		private final long key2;
		@Override
		public int hashCode() {
			long xor = key1 ^ key2;
			return (int)(xor ^ (xor >>> 32));
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Hash128bit other = (Hash128bit) obj;
			if (key1 != other.key1)
				return false;
			if (key2 != other.key2)
				return false;
			return true;
		}
		private static final HashFunction hf = Hashing.murmur3_128();
		public Hash128bit(String key) {
			HashCode hc = hf.newHasher().putString(key, StandardCharsets.UTF_8).hash();
			ByteBuffer bb = ByteBuffer.wrap(hc.asBytes());
			key1 = bb.getLong();
			key2 = bb.getLong();
		}
		@Override
		public String toString() {
			return Long.toHexString(key1) + Long.toHexString(key2); 
		}
	}
}
