package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import java.util.List;
import java.util.NavigableMap;

import com.google.common.collect.Maps;

/**
 * Calculates the coordinate of a genomic position based on the concatenation
 * of all chromosomes in sequence dictionary order
 * @author Daniel Cameron
 *
 */
public class LinearGenomicCoordinate {
	private final SAMSequenceDictionary dictionary;
	private long[] offset;
	private NavigableMap<Long, Integer> referenceIndices;
	/**
	 * Create a coordinate lookup from the given dictionary
	 * @param dictionary sequence dictionary
	 * @param padding number of bases to pad at start and end of each chromosome
	 */
	public LinearGenomicCoordinate(SAMSequenceDictionary dictionary, long padding) {
		this.dictionary = dictionary;
		generateLookups(dictionary.getSequences(), padding);
	}
	public LinearGenomicCoordinate(SAMSequenceDictionary dictionary) {
		this(dictionary, 0);
	}
	private void generateLookups(List<SAMSequenceRecord> sequences, long padding) {
		offset = new long[sequences.size()];
		long cumsum = 0;
		for (int i = 0; i < sequences.size(); i++) {
			cumsum += padding;
			offset[i] = cumsum;
			cumsum += sequences.get(i).getSequenceLength() + padding;
			if (sequences.get(i).getSequenceLength() == 0) {
				throw new IllegalArgumentException(String.format("Missing length for contig %s", sequences.get(i).getSequenceName()));
			}
		}
		referenceIndices = Maps.newTreeMap();
		for (int i = 0; i < sequences.size(); i++) {
			referenceIndices.put(offset[i] - padding, i);
		}
		// add sentinels so lookup always succeeds
		referenceIndices.put(Long.MIN_VALUE, -1);
		referenceIndices.put(cumsum + 1, -1);
	}
	public long getLinearCoordinate(int referenceIndex, long pos) {
		if (referenceIndex < 0 || referenceIndex >= offset.length) {
			return -1;
		}
		return offset[referenceIndex] + pos;
	}
	public long getLinearCoordinate(String chr, long pos) {
		return getLinearCoordinate(dictionary.getSequenceIndex(chr), pos);
	}
	public long getStartLinearCoordinate(BreakendSummary bp) {
		return getLinearCoordinate(bp.referenceIndex, bp.start);
	}
	public long getStartLinearCoordinate(SAMRecord r) {
		return getLinearCoordinate(r.getReferenceIndex(), r.getAlignmentStart());
	}
	public long getEndLinearCoordinate(BreakendSummary bp) {
		return getLinearCoordinate(bp.referenceIndex, bp.end);
	}
	public long getEndLinearCoordinate(SAMRecord r) {
		return getLinearCoordinate(r.getReferenceIndex(), r.getAlignmentEnd());
	}
	public int getReferenceIndex(long linearCoordinate) {
		return referenceIndices.floorEntry(linearCoordinate - 1).getValue();
	}
	public int getReferencePosition(long linearCoordinate) {
		int referenceIndex = getReferenceIndex(linearCoordinate);
		if (referenceIndex < 0) return -1;
		return (int)(linearCoordinate - getLinearCoordinate(referenceIndex, 0));
	}
}
