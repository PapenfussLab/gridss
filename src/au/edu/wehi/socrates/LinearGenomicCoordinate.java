package au.edu.wehi.socrates;

import java.util.List;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Calculates the coordinate of a genomic position based on the concatenation
 * of all chromosomes in sequence dictionary order
 * @author Daniel Cameron
 *
 */
public class LinearGenomicCoordinate {
	private SAMSequenceDictionary dictionary;
	private long[] offset;
	/**
	 * Create a coordinate lookup from the given dictionary
	 * @param dictionary sequence dictionary
	 * @param number of bases to pad around chromosomes 
	 */
	public LinearGenomicCoordinate(SAMSequenceDictionary dictionary, int buffer) {
		this.dictionary = dictionary;
		this.offset = generateLookup(dictionary.getSequences(), buffer);
	}
	public LinearGenomicCoordinate(SAMSequenceDictionary dictionary) {
		this(dictionary, 0);
	}
	private static long[] generateLookup(List<SAMSequenceRecord> sequences, int buffer) {
		long[] lookup = new long[sequences.size()];
		long cumsum = buffer;
		for (int i = 0; i < sequences.size(); i++) {
			lookup[i] = cumsum;
			cumsum += sequences.get(i).getSequenceLength() + buffer;
			if (sequences.get(i).getSequenceLength() == 0) {
				throw new IllegalArgumentException(String.format("Missing length for contig %s", sequences.get(i).getSequenceName()));
			}
		}
		return lookup;
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
	public long getStartLinearCoordinate(BreakpointLocation bp) {
		return getLinearCoordinate(bp.referenceIndex, bp.start);
	}
}
