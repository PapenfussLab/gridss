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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;

public class GreedyAllocationCache {
	private static final Comparator<ChimericAlignment> ByPosition = Comparator.<ChimericAlignment, Integer>comparing(ca -> ca.pos)
				.thenComparing(ca -> ca.isNegativeStrand)
				.thenComparing(ca -> ca.rname)
				.thenComparing(ca -> ca.cigar.toString());
	/**
	 * Gets a string representing the full read alignment including split alignments
	 * @param r record
	 * @return read alignment string
	 */
	protected String getReadAlignment(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return "";
		List<ChimericAlignment> ca = Lists.newArrayList(ChimericAlignment.getChimericAlignments(r));
		ca.add(new ChimericAlignment(r));
		return getAlignmentString(ca);
	}
	/**
	 * Gets a string representing the primary alignment of the read pair
	 * @param r
	 * @return
	 */
	protected Hash128bit getReadPairAlignment(SAMRecord r) {
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
	protected static class Node {
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
	protected void put(HashMap<Hash128bit, Node> bestAlignment, Hash128bit readName, Hash128bit alignment, Hash128bit event, float score) {
		if (bestAlignment == null) return;
		Node existing = bestAlignment.get(readName);
		if (existing == null || existing.score < score) {
			bestAlignment.put(readName, new Node(event, score, alignment));
		}
	}

	protected boolean isBestAlignment(HashMap<Hash128bit, Node> bestAlignment, Hash128bit readName, Hash128bit alignment) {
		if (bestAlignment == null) return true;
		Node node = bestAlignment.get(readName);
		return node != null && Objects.equals(alignment, node.alignment);
	}
	protected static class Hash128bit {
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
		public boolean equals(Hash128bit obj) {
			return obj != null && key1 == obj.key1 && key2 == obj.key2;
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