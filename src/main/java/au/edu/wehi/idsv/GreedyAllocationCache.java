package au.edu.wehi.idsv;

import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
	protected static final int INITIAL_HASH_MAP_SIZE = 1 << 20;
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
	protected static class AlignmentScoreNode extends Hash128bit {
		public AlignmentScoreNode(float score, Hash128bit alignment) {
			super(alignment.key1, alignment.key2);
			this.score = score;
		}
		private final float score;
		public Hash128bit getAlignment() { return this; } 
		public float getScore() { return score; }
		public String toString() {
			return String.format("(%f,%s)", getAlignment(), getScore());
		}
	}
	protected static class EventAlignmentScoreNode extends AlignmentScoreNode {
		public EventAlignmentScoreNode(Hash128bit event, float score, Hash128bit alignment) {
			super(score, alignment);
			this.eventKey1 = event.key1;
			this.eventKey2 = event.key2;
		}
		private final long eventKey1;
		private final long eventKey2;
		public Hash128bit getEvent() { return new Hash128bit(eventKey1, eventKey2); }
		public String toString() {
			return String.format("(%s,%f,%s)", getEvent(), getScore(), getAlignment());
		}
	}
	protected static class EventScoreNode extends Hash128bit {
		public EventScoreNode(Hash128bit event, float score) {
			super(event.key1, event.key2);
			this.score = score;
		}
		private final float score;
		public Hash128bit getEvent() { return this; } 
		public float getScore() { return score; }
		public String toString() {
			return String.format("(,%s,%f)", getEvent(), getScore());
		}
	}
	protected void put(Map<Hash128bit, EventAlignmentScoreNode> bestAlignment, Hash128bit readName, Hash128bit alignment, Hash128bit event, float score) {
		if (bestAlignment == null) return;
		EventAlignmentScoreNode notInLookup = new EventAlignmentScoreNode(event, score, alignment);
		EventAlignmentScoreNode inLookup = bestAlignment.get(readName);
		// this check is a loop as when there are multiple threads writing to the same key
		// the score that we wrote could be worse since the get and put are two separate (atomic)
		// operations.
		while (notInLookup != null && (inLookup == null || inLookup.getScore() < notInLookup.getScore())) {
			EventAlignmentScoreNode toPut = notInLookup;
			inLookup = toPut;
			notInLookup = bestAlignment.put(readName, toPut);
		}
	}
	protected void put(Map<Hash128bit, AlignmentScoreNode> bestAlignment, Hash128bit readName, Hash128bit alignment, float score) {
		// (duplicates EventAlignmentScoreNode logic)
		if (bestAlignment == null) return;
		AlignmentScoreNode notInLookup = new AlignmentScoreNode(score, alignment);
		AlignmentScoreNode inLookup = bestAlignment.get(readName);
		while (notInLookup != null && (inLookup == null || inLookup.getScore() < notInLookup.getScore())) {
			AlignmentScoreNode toPut = notInLookup;
			inLookup = toPut;
			notInLookup = bestAlignment.put(readName, toPut);
		}
	}
	protected void put(HashMap<Hash128bit, EventScoreNode> bestEvent, Hash128bit readName, Hash128bit event, float score) {
		if (bestEvent == null) return;
		EventScoreNode notInLookup = new EventScoreNode(event, score);
		EventScoreNode inLookup = bestEvent.get(readName);
		while (notInLookup != null && (inLookup == null || inLookup.getScore() < notInLookup.getScore())) {
			EventScoreNode toPut = notInLookup;
			inLookup = toPut;
			notInLookup = bestEvent.put(readName, toPut);
		}
	}
	protected boolean isBestAlignment(Map<Hash128bit, EventAlignmentScoreNode> bestAlignment, Hash128bit readName, Hash128bit alignment) {
		if (bestAlignment == null) return true;
		EventAlignmentScoreNode node = bestAlignment.get(readName);
		return node != null && Objects.equals(alignment, node.getAlignment());
	}
	protected static class Hash128bit {
		protected final long key1;
		protected final long key2;
		@Override
		public int hashCode() {
			long xor = key1 ^ key2;
			return (int)(xor ^ (xor >>> 32));
		}
		@Override
		public boolean equals(Object obj) {
			if (obj instanceof Hash128bit) {
				return equals((Hash128bit)obj);
			}
			return false;
		}
		public boolean equals(Hash128bit obj) {
			return obj != null && key1 == obj.key1 && key2 == obj.key2;
		}
		private static final HashFunction hf = Hashing.murmur3_128();
		public Hash128bit(String key) {
			HashCode hc = hf.newHasher().putString(key, StandardCharsets.UTF_8).hash();
			ByteBuffer bb = ByteBuffer.wrap(hc.asBytes());
			this.key1 = bb.getLong();
			this.key2 = bb.getLong();
		}
		public Hash128bit(long key1, long key2) {
			this.key1 = key1;
			this.key2 = key2;
		}
		@Override
		public String toString() {
			return Long.toHexString(key1) + Long.toHexString(key2); 
		}
	}
}