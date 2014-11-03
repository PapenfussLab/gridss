package au.edu.wehi.idsv.util;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;

import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;

/**
 * Buffers entire reference to enable efficient random lookup of sequences
 * @author cameron.d
 *
 */
public class BufferedReferenceSequenceFile implements ReferenceSequenceFile {
	private static final Log log = Log.getInstance(BufferedReferenceSequenceFile.class);
	private final ReferenceSequenceFile underlying;
	private final Map<String, ReferenceSequence> cache = Maps.newHashMap();
	public BufferedReferenceSequenceFile(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
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
	@Override
	public ReferenceSequence getSequence(String contig) {
		if (!cache.containsKey(contig)) {
			cache.put(contig, underlying.getSequence(contig));
			log.debug("Cached reference genome contig ", contig);
		}
		return cache.get(contig);
	}
	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
        int length = (int)(stop - start + 1);
		ReferenceSequence fullContig = getSequence(contig);
		if (length > fullContig.length()) {
			throw new IllegalArgumentException("subsequence out of contig bounds");
		}
		if (start > stop + 1) {
			throw new IllegalArgumentException("start after stop");
		}
		byte[] target = new byte[length];
		System.arraycopy(fullContig.getBases(), (int) (start - 1), target, 0, target.length);
		return new ReferenceSequence(fullContig.getName(), fullContig.getContigIndex(), target);
	}
	@Override
	public void close() throws IOException {
		underlying.close();	
	}
}
