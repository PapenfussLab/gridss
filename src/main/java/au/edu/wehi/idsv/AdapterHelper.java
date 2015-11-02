package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.List;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedKmerList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Provides adapter related utilities
 * @author cameron.d
 *
 */
public class AdapterHelper {
	private static final int MAX_ADAPTER_SOFT_CLIP_MICROHOMOLOGY_LENGTH = 6;
	private List<String> adapterSequences;
	private int k;
	private long[] kmers;
	public AdapterHelper(String[] adapters) {
		if (adapters == null) adapters = new String[0];
		this.adapterSequences = ImmutableList.copyOf(adapters);
		if (hasAdapters()) {
			// cache stuff
			k = adapterSequences.get(0).length();
			kmers = new long[adapterSequences.size() * 2];
			for (int i = 0; i < adapterSequences.size(); i++) {
				kmers[2 * i] = KmerEncodingHelper.picardBaseToEncoded(k, adapterSequences.get(i).getBytes(StandardCharsets.US_ASCII));
				kmers[2 * i + 1] = KmerEncodingHelper.reverseComplement(k, kmers[2 * i]);
			}
		}
	}
	private boolean hasAdapters() { return adapterSequences.size() > 0; }
	public String[] getAdapterSequences() {
		return adapterSequences.toArray(new String[0]);
	}
	public boolean containsAdapter(SAMRecord record) {
		if (!hasAdapters()) return false;
		if (record.getReadLength() < k) return false;
		PackedKmerList list = new PackedKmerList(k, record.getReadBases(), null, false, false);
		for (int i = 0; i < list.length(); i++) {
			long kmer = list.kmer(i);
			for (int j = 0; j < kmers.length; j++) {
				if (kmer == kmers[j]) {
					return true;
				}
			}
		}
		return false;
	}
	/**
	 * Determine whether this soft clip is cause by read-through into adapter sequence 
	 * @param p soft clip parameter
	 * @return true, if the soft clip is due to adapter sequence, false otherwise
	 */
	public boolean isAdapterSoftClip(SoftClipEvidence e) {
		if (!hasAdapters()) return false;
		BreakendDirection direction = e.getBreakendSummary().direction;
		SAMRecord record = e.getSAMRecord();
		if (direction == BreakendDirection.Forward && record.getReadNegativeStrandFlag()) return false;
		if (direction == BreakendDirection.Backward && !record.getReadNegativeStrandFlag()) return false;
		for (String adapter : adapterSequences) {
			if (scMatchesAdapterFR(record, direction, e.getSoftClipLength(), adapter)) {
				return true;
			}
		}
		return false;
	}
	/**
	 * Checks the soft clip against the given adapter
	 * @param adapter
	 * @return
	 */
	private boolean scMatchesAdapterFR(SAMRecord record, BreakendDirection direction, int scLen, String adapter) {
		if (direction == BreakendDirection.Forward) {
			// match soft clip
			for (int i = 0; i <= MAX_ADAPTER_SOFT_CLIP_MICROHOMOLOGY_LENGTH; i++) {
				if (matchesAdapterSequence(adapter, record.getReadBases(), record.getReadLength() - scLen - i, 1, false)) {
					return true;
				}
			}
		} else {
			for (int i = 0; i <= MAX_ADAPTER_SOFT_CLIP_MICROHOMOLOGY_LENGTH; i++) {
				if (matchesAdapterSequence(adapter, record.getReadBases(), scLen + i - 1, -1, true)) {
					return true;
				}
			}
		}
		return false;
	}
	private static boolean matchesAdapterSequence(String adapter, byte[] read, int readStartOffset, int readDirection, boolean complementAdapter) {
		boolean nonzeroSize = false;
		for (int i = 0; readStartOffset + i * readDirection < read.length && readStartOffset + i * readDirection  >= 0 && i < adapter.length(); i++) {
			byte readBase = read[readStartOffset + i * readDirection];
			byte adapterBase = (byte)adapter.charAt(i);
			if (complementAdapter) adapterBase = SequenceUtil.complement(adapterBase);
			if (SequenceUtil.isValidBase(readBase) && readBase != adapterBase) {
				return false;
			}
			nonzeroSize = true;
		}
		return nonzeroSize; // no bases compared = no adapter match 
	}
}
