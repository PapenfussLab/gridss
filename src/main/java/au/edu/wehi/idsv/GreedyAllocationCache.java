package au.edu.wehi.idsv;

import java.io.Closeable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.ignite.Ignite;
import org.apache.ignite.IgniteCache;
import org.apache.ignite.IgniteLogger;
import org.apache.ignite.Ignition;
import org.apache.ignite.binary.BinaryObjectException;
import org.apache.ignite.binary.BinaryRawReader;
import org.apache.ignite.binary.BinaryRawWriter;
import org.apache.ignite.binary.BinaryReader;
import org.apache.ignite.binary.BinaryWriter;
import org.apache.ignite.binary.Binarylizable;
import org.apache.ignite.cache.CacheMemoryMode;
import org.apache.ignite.configuration.CacheConfiguration;
import org.apache.ignite.configuration.IgniteConfiguration;

import com.google.common.collect.Lists;
import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Log;
import net.openhft.chronicle.bytes.BytesIn;
import net.openhft.chronicle.bytes.BytesMarshallable;
import net.openhft.chronicle.bytes.BytesOut;
import net.openhft.chronicle.map.ChronicleMap;

public abstract class GreedyAllocationCache implements Closeable {
	//private static final Log log = Log.getInstance(GreedyAllocationCache.class);
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
	protected void putEventAlignmentScoreNode(GreedyAllocationCacheLookup<EventAlignmentScoreNode> bestAlignment, Hash96bit readName, Hash96bit alignment, Hash96bit event, float score) {
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
	protected void putAlignmentScoreNode(GreedyAllocationCacheLookup<AlignmentScoreNode> bestAlignment, Hash96bit readName, Hash96bit alignment, float score) {
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
	protected void putEventScoreNode(GreedyAllocationCacheLookup<EventScoreNode> bestEvent, Hash96bit readName, Hash96bit event, float score) {
		if (bestEvent == null) return;
		EventScoreNode notInLookup = new EventScoreNode(event, score);
		EventScoreNode inLookup = bestEvent.get(readName);
		while (notInLookup != null && (inLookup == null || inLookup.getScore() < notInLookup.getScore())) {
			EventScoreNode toPut = notInLookup;
			inLookup = toPut;
			notInLookup = bestEvent.put(readName, toPut);
		}
	}
	protected boolean isBestAlignment(GreedyAllocationCacheLookup<EventAlignmentScoreNode> bestAlignment, Hash96bit readName, Hash96bit alignment) {
		if (bestAlignment == null) return true;
		EventAlignmentScoreNode node = bestAlignment.get(readName);
		return node != null && Objects.equals(alignment, node.getAlignment());
	}
	/**
	 * Gets a string representing the primary alignment of the read pair
	 * @param r
	 * @return
	 */
	protected Hash96bit getReadPairAlignment(SAMRecord r) {
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
		return new Hash96bit(getAlignmentString(ca));
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
	protected static interface GreedyAllocationCacheLookup<T extends Hash96bit> extends Closeable {
		T get(Hash96bit key);
		T put(Hash96bit key, T value);
	}
	protected static <T extends Hash96bit> GreedyAllocationCacheLookup<T> createLookup(String name, Class<T> clazz, long size) {
		//return new HashMapLookup<T>();
		//return new ApacheIgniteLookup<T>(name);
		return new OpenHFTLookup<T>(name, clazz, size);
	}
	/**
	 * Java HashMap lookup. This stores all values on the java heap which requireds
	 * results a gigantic (>400GB for WGS) heap
	 */
	protected static class HashMapLookup<T extends Hash96bit> implements GreedyAllocationCacheLookup<T> {
		private static final int INITIAL_LOOKUP_SIZE = 1 << 20;
		private Map<Hash96bit, T> map;
		public HashMapLookup() {
			//map = new HashMap<>(INITIAL_LOOKUP_SIZE);
			map = new ConcurrentHashMap<Hash96bit, T>(INITIAL_LOOKUP_SIZE, 0.75f, 64);
		}
		@Override
		public T get(Hash96bit key) {
			return map.get(key);
		}
		@Override
		public T put(Hash96bit key, T value) {
			return map.put(key, value);
		}
		@Override
		public void close() throws IOException {
			map = null;
		}
	}
	protected static class ApacheIgniteLookup<T extends Hash96bit> implements GreedyAllocationCacheLookup<T> {
		private static final Log log = Log.getInstance(ApacheIgniteLookup.class);
		static {
			IgniteConfiguration cfg = new IgniteConfiguration();
			cfg.setGridName("GreedyAllocationCache");
			cfg.setGridLogger(new Logger());
			ignite = Ignition.start(cfg);
		}
		private static Ignite ignite;
		private String name;
		private IgniteCache<Hash96bit, T> cache;
		public ApacheIgniteLookup(String name) {
			CacheConfiguration<Hash96bit, T> cacheCfg = new CacheConfiguration<>();
			cacheCfg.setMemoryMode(CacheMemoryMode.OFFHEAP_TIERED); // store everything off-heap
			cacheCfg.setOffHeapMaxMemory(0); // Set off-heap memory unlimited
			//cacheCfg.setEvictionPolicy(new FifoEvictionPolicy<>(INITIAL_LOOKUP_SIZE, INITIAL_LOOKUP_SIZE / 10));
			cacheCfg.setName(name);
			this.name = name;
			this.cache = ignite.getOrCreateCache(cacheCfg);
		}
		@Override
		public T get(Hash96bit key) {
			return cache.get(key);
		}
		@Override
		public T put(Hash96bit key, T value) {
			return cache.getAndPut(key, value);
		}
		@Override
		public void close() throws IOException {
			cache.destroy();
			ignite.destroyCache(name);
		}
		private static class Logger implements IgniteLogger {
			public IgniteLogger getLogger(Object ctgr) { return this; }
			public void trace(String msg) { }
			public void debug(String msg) { }
			public void info(String msg) { log.info(msg); }
			public void warning(String msg) { log.warn(msg); }
			public void warning(String msg, Throwable e) { log.warn(e, msg); }
			public void error(String msg) { log.error(msg); }
			public void error(String msg, Throwable e) { log.error(e, msg); }
			public boolean isTraceEnabled() { return false; }
			public boolean isDebugEnabled() { return false; }
			public boolean isInfoEnabled() { return true; }
			public boolean isQuiet() { return true; }
			public String fileName() { return null; }
		}
	}
	protected static class OpenHFTLookup<T extends Hash96bit> implements GreedyAllocationCacheLookup<T> {
		private ChronicleMap<Hash96bit, T> map;
		private ChronicleMap<Hash96bit, T> create(String name, Class<T> clazz, long size) {
			try {
				return ChronicleMap
					    .of(Hash96bit.class, clazz)
					    .name(name)
					    .constantKeySizeBySample(new Hash96bit())
					    .constantValueSizeBySample(clazz.newInstance())
					    .entries(size)
					    .create();
			} catch (InstantiationException | IllegalAccessException e) {
				throw new RuntimeException(e);
			}
		}
		public OpenHFTLookup(String name, Class<T> clazz, long size) {
			this.map = create(name, clazz, size);
		}
		@Override
		public void close() throws IOException {
			map.close();
		}
		@Override
		public T get(Hash96bit key) {
			return map.get(key);
		}
		@Override
		public T put(Hash96bit key, T value) {
			return map.put(key, value);
		}
	}
	protected static class Hash96bit implements Binarylizable, BytesMarshallable {
		private static final HashFunction hf = Hashing.murmur3_128();
		protected long key1;
		protected int key2;
		public Hash96bit() { }
		public Hash96bit(String key) {
			HashCode hc = hf.newHasher().putString(key, StandardCharsets.UTF_8).hash();
			ByteBuffer bb = ByteBuffer.wrap(hc.asBytes());
			this.key1 = bb.getLong();
			this.key2 = bb.getInt();
		}
		public Hash96bit(Hash96bit hash) {
			this.key1 = hash.key1;
			this.key2 = hash.key2;
		}
		public Hash96bit(long key1, int key2) {
			this.key1 = key1;
			this.key2 = key2;
		}
		@Override
		public int hashCode() {
			long xor = key1 ^ key2;
			return (int)(xor ^ (xor >>> 32));
		}
		@Override
		public boolean equals(Object obj) {
			if (obj instanceof Hash96bit) {
				return equals((Hash96bit)obj);
			}
			return false;
		}
		public boolean equals(Hash96bit obj) {
			return obj != null && key1 == obj.key1 && key2 == obj.key2;
		}
		@Override
		public String toString() {
			return Long.toHexString(key1) + Long.toHexString(key2); 
		}
		@Override
		public void	readBinary(BinaryReader reader) {
			BinaryRawReader brr = reader.rawReader();
			key1 = brr.readLong();
			key2 = brr.readInt();
		}
		@Override
		public void writeBinary(BinaryWriter writer) throws BinaryObjectException {
			BinaryRawWriter brw = writer.rawWriter();
			brw.writeLong(key1);
			brw.writeInt(key2);
		}
		@SuppressWarnings("rawtypes")
		@Override
		public void readMarshallable(BytesIn bytes) {
			key1 = bytes.readLong();
			key2 = bytes.readInt();
		}
		@SuppressWarnings("rawtypes")
		@Override
		public void writeMarshallable(BytesOut bytes) {
			bytes.writeLong(key1);
			bytes.writeInt(key2);
		}
	}
	protected static abstract class HashScoreNode extends Hash96bit {
		public HashScoreNode() { }
		public HashScoreNode(Hash96bit hash, float score) {
			super(hash);
			this.score = score;
		}
		protected float score; 
		public float getScore() { return score; }
		@Override
		public void readBinary(BinaryReader reader) {
			BinaryRawReader brr = reader.rawReader();
			key1 = brr.readLong();
			key2 = brr.readInt();
			score = brr.readFloat();
		}
		@Override
		public void writeBinary(BinaryWriter writer) throws BinaryObjectException {
			BinaryRawWriter brw = writer.rawWriter();
			brw.writeLong(key1);
			brw.writeInt(key2);
			brw.writeFloat(score);
		}
		@SuppressWarnings("rawtypes")
		@Override
		public void readMarshallable(BytesIn bytes) {
			key1 = bytes.readLong();
			key2 = bytes.readInt();
			score = bytes.readFloat();
		}
		@SuppressWarnings("rawtypes")
		@Override
		public void writeMarshallable(BytesOut bytes) {
			bytes.writeLong(key1);
			bytes.writeInt(key2);
			bytes.writeFloat(score);
		}
	}
	protected static class AlignmentScoreNode extends HashScoreNode {
		public AlignmentScoreNode() { }
		public AlignmentScoreNode(float score, Hash96bit alignment) {
			super(alignment, score);
		}
		public Hash96bit getAlignment() { return this; } 
		public String toString() {
			return String.format("(%f,%s)", getAlignment(), getScore());
		}
	}
	protected static class EventScoreNode extends HashScoreNode {
		public EventScoreNode() { }
		public EventScoreNode(Hash96bit event, float score) {
			super(event, score);
		}
		public Hash96bit getEvent() { return this; } 
		public String toString() {
			return String.format("(,%s,%f)", getEvent(), getScore());
		}
	}
	protected static class EventAlignmentScoreNode extends AlignmentScoreNode {
		public EventAlignmentScoreNode() { }
		public EventAlignmentScoreNode(Hash96bit event, float score, Hash96bit alignment) {
			super(score, alignment);
			this.eventKey1 = event.key1;
			this.eventKey2 = event.key2;
		}
		protected long eventKey1;
		protected int eventKey2;
		public Hash96bit getEvent() { return new Hash96bit(eventKey1, eventKey2); }
		public String toString() {
			return String.format("(%s,%f,%s)", getEvent(), getScore(), getAlignment());
		}
		@Override
		public void readBinary(BinaryReader reader) {
			BinaryRawReader brr = reader.rawReader();
			key1 = brr.readLong();
			key2 = brr.readInt();
			eventKey1 = brr.readLong();
			eventKey2 = brr.readInt();
			score = brr.readFloat();
		}
		@Override
		public void writeBinary(BinaryWriter writer) throws BinaryObjectException {
			BinaryRawWriter brw = writer.rawWriter();
			brw.writeLong(key1);
			brw.writeInt(key2);
			brw.writeLong(eventKey1);
			brw.writeInt(eventKey2);
			brw.writeFloat(score);
		}
		@SuppressWarnings("rawtypes")
		@Override
		public void readMarshallable(BytesIn bytes) {
			key1 = bytes.readLong();
			key2 = bytes.readInt();
			eventKey1 = bytes.readLong();
			eventKey2 = bytes.readInt();
			score = bytes.readFloat();
		}
		@SuppressWarnings("rawtypes")
		@Override
		public void writeMarshallable(BytesOut bytes) {
			bytes.writeLong(key1);
			bytes.writeInt(key2);
			bytes.writeLong(eventKey1);
			bytes.writeInt(eventKey2);
			bytes.writeFloat(score);
		}
	}
}