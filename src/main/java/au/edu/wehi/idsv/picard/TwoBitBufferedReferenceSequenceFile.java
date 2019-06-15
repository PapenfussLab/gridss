package au.edu.wehi.idsv.picard;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedSequence;
import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;

import java.io.*;
import java.nio.MappedByteBuffer;
import java.nio.file.Files;
import java.time.OffsetTime;
import java.util.BitSet;
import java.util.stream.Stream;

/**
 * 2bit encodes and buffers the entire reference to enable efficient random lookup of small subsequences
 * @author Daniel Cameron
 *
 */
public class TwoBitBufferedReferenceSequenceFile implements ReferenceSequenceFile, ReferenceLookup {
	private static final Log log = Log.getInstance(TwoBitBufferedReferenceSequenceFile.class);
	private final ReferenceSequenceFile underlying;
	private final PackedReferenceSequence[] referenceIndexLookup;
	private final File cacheFile;
	private boolean attemptedSave = false;
	/**
	 * Cached contigs
	 */
	private ImmutableMap<String, PackedReferenceSequence> cache = ImmutableMap.of();
	public TwoBitBufferedReferenceSequenceFile(ReferenceSequenceFile underlying) {
		this(underlying, null);
	}
	public TwoBitBufferedReferenceSequenceFile(ReferenceSequenceFile underlying, File cache) {
		this.underlying = underlying;
		this.referenceIndexLookup = new PackedReferenceSequence[underlying.getSequenceDictionary().getSequences().size()];
		this.cacheFile = cache;
	}
	public byte getBase(int referenceIndex, int position) {
		PackedReferenceSequence seq = referenceIndexLookup[referenceIndex];
		if (seq == null) {
			seq = addToCache(underlying.getSequenceDictionary().getSequence(referenceIndex).getSequenceName());
		}
		if (seq.ambiguous.get(position - 1)) {
			return 'N';
		}
		return seq.get(position - 1);
	}
	public synchronized void load(File file) {
		ImmutableMap.Builder<String, PackedReferenceSequence> builder = ImmutableMap.<String, PackedReferenceSequence>builder();
		try(FileInputStream fis = new FileInputStream(file)) {
			try(ObjectInputStream ois = new ObjectInputStream(fis)) {
				for (int i =  0; i < referenceIndexLookup.length; i++) {
					referenceIndexLookup[i] = (PackedReferenceSequence)ois.readObject();
					builder.put(referenceIndexLookup[i].name, referenceIndexLookup[i]);
				}
			}
			cache = builder.build();
		} catch (Exception e) {
			log.error("Error loading reference genome from cache " + file, e);
		}
	}
	public synchronized void save(File file) {
		if (file.exists()) {
			throw new IllegalArgumentException(file + " already exists");
		}
		// Ensure the lookup is fully populated
		underlying.getSequenceDictionary()
				.getSequences()
				.stream()
				.map(s -> s.getSequenceName())
				.forEach(s -> cacheLoad(s));
		try(FileOutputStream fos = new FileOutputStream(file)) {
			try(ObjectOutputStream oos = new ObjectOutputStream(fos)) {
				for (int i =  0; i < referenceIndexLookup.length; i++) {
					oos.writeObject(referenceIndexLookup[i]);
				}
			}
		} catch (Exception e) {
			log.error("Error saving reference genome to cache file " + file, e);
			try {
				if (file.exists()) {
					Files.delete(file.toPath());
				}
			} catch (IOException e1) {
				// swallow recovery exception
			}
		}
	}
	private static class PackedReferenceSequence extends PackedSequence implements Serializable {
		private final String name;
	    private final int contigIndex;
	    private final long length;
	    private final BitSet ambiguous;
		public PackedReferenceSequence(ReferenceSequence seq) {
			super(seq.getBases(), false, false);
			this.name = seq.getName();
			this.contigIndex = seq.getContigIndex();
			this.length = seq.length();
			this.ambiguous = new BitSet(seq.length());
			byte[] seqBases = seq.getBases();
			for (int i = 0; i < length; i++) {
				if (KmerEncodingHelper.isAmbiguous(seqBases[i])) {
					ambiguous.set(i);
				}
			}
		}
		public ReferenceSequence getSequence() {
			return getSubsequenceAt(1, length);
		}
		public ReferenceSequence getSubsequenceAt(long start, long stop) {
			int length = (int)(stop - start + 1);
			ReferenceSequence seq = new ReferenceSequence(name, contigIndex, getBytes((int)(start - 1), length));
			byte[] seqBases = seq.getBases();
			for (int i = 0; i < length; i++) {
				if (ambiguous.get((int)start - 1 + i)) {
					seqBases[i] = 'N';
				}
			}
			return seq;
		}
	}
	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return underlying.getSequenceDictionary();
	}
	@Override
	public ReferenceSequence nextSequence() {
		return underlying.nextSequence();
	}
	@Override
	public void reset() {
		underlying.reset();
	}
	@Override
	public boolean isIndexed() {
		return underlying.isIndexed();
	}
	/**
	 * Updates the cache to include the new contig
	 * @param contig
	 */
	private synchronized PackedReferenceSequence addToCache(String contig) {
		if (cacheFile != null) {
			if (cacheFile.exists()) {
				log.info("Loading reference genome from cache " + cacheFile);
				load(cacheFile);
				log.info("Loading reference genome complete");
			} else if (!attemptedSave) {
				attemptedSave = true;
				if (!cacheFile.getParentFile().canWrite()) {
					log.warn("Cannot write to " + cacheFile + " not persisting 2bit compressed reference genome cache");
				} else {
					log.info("Saving reference genome cache to " + cacheFile);
					save(cacheFile);
					log.info("Saving reference genome cache complete");
				}
			}
		}
		return cacheLoad(contig);
	}
	private synchronized PackedReferenceSequence cacheLoad(String contig) {
		PackedReferenceSequence seq = cache.get(contig);
		if (seq != null) {
			// already populated by another thread while we were waiting to enter
			// this synchronized block
			return seq;
		}
		log.info("Caching reference genome contig " + contig);
		ReferenceSequence fullContigSequence = underlying.getSequence(contig);
		seq = new PackedReferenceSequence(fullContigSequence);
		cache = ImmutableMap.<String, PackedReferenceSequence>builder()
				.putAll(cache)
				.put(contig, seq)
				.build();
		referenceIndexLookup[underlying.getSequenceDictionary().getSequence(contig).getSequenceIndex()] = seq;
		return seq;
	}
	@Override
	public ReferenceSequence getSequence(String contig) {
		PackedReferenceSequence seq = cache.get(contig);
		if (seq == null) {
			seq = addToCache(contig);
		}
		return seq.getSequence();
	}
	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
        PackedReferenceSequence seq = cache.get(contig);
		if (seq == null) {
			seq = addToCache(contig);
		}
		return seq.getSubsequenceAt(start, stop);
	}
	@Override
	public void close() throws IOException {
		underlying.close();	
	}
}
