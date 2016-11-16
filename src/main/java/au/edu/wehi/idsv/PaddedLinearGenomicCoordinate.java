package au.edu.wehi.idsv;

import java.util.List;
import java.util.NavigableMap;

import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Calculates the coordinate of a genomic position based on the concatenation
 * of all chromosomes in sequence dictionary order with padding between
 * chromosomes
 * @author Daniel Cameron
 *
 */
public class PaddedLinearGenomicCoordinate implements LinearGenomicCoordinate {
	private final SAMSequenceDictionary dictionary;
	private long[] offset;
	private NavigableMap<Long, Integer> referenceIndices;
	/**
	 * Create a coordinate lookup from the given dictionary
	 * @param dictionary sequence dictionary
	 * @param padding number of bases to pad at start and end of each chromosome
	 * @param padding between chromosomes include sequence bases such at each chromosome starts at an integer multiple of padding
	 */
	public PaddedLinearGenomicCoordinate(SAMSequenceDictionary dictionary, long padding, boolean fixedPadding) {
		this.dictionary = dictionary;
		generateLookups(dictionary.getSequences(), padding, fixedPadding);
	}
	/**
	 * Create a coordinate lookup from the given dictionary
	 * @param dictionary sequence dictionary
	 * @param padding number of bases to pad at start and end of each chromosome
	 */
	public PaddedLinearGenomicCoordinate(SAMSequenceDictionary dictionary, long padding) {
		this(dictionary, padding, false);
	}
	public PaddedLinearGenomicCoordinate(SAMSequenceDictionary dictionary) {
		this(dictionary, 0);
	}
	private void generateLookups(List<SAMSequenceRecord> sequences, long padding, boolean fixedPadding) {
		if (fixedPadding) {
			generateFixedLookups(sequences, padding);
		} else {
			generatePaddedLookups(sequences, padding);
		}
	}
	private void generatePaddedLookups(List<SAMSequenceRecord> sequences, long padding) {
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
		referenceIndices.put(cumsum, -1);
	}
	private void generateFixedLookups(List<SAMSequenceRecord> sequences, long width) {
		offset = new long[sequences.size()];
		for (int i = 0; i < sequences.size(); i++) {
			if (sequences.get(i).getSequenceLength() > width) {
				throw new IllegalArgumentException(String.format("Length %dbp of %s  is longer than fixed padding size of %d", sequences.get(i).getSequenceLength(), sequences.get(i).getSequenceName(), width));
			}
			if (sequences.get(i).getSequenceLength() == 0) {
				throw new IllegalArgumentException(String.format("Missing length for contig %s", sequences.get(i).getSequenceName()));
			}
			offset[i] = (i + 1) * width;
		}
		referenceIndices = Maps.newTreeMap();
		for (int i = 0; i < sequences.size(); i++) {
			referenceIndices.put(offset[i] - width / 2, i);
		}
		referenceIndices.put(0L, 0); // first gets a full width of starting padding
		// add sentinels so lookup always succeeds
		referenceIndices.put(Long.MIN_VALUE, -1);
		referenceIndices.put((sequences.size() + 2) * width - 1, -1); // end gets a full width of ending padding
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getLinearCoordinate(int, long)
	 */
	@Override
	public long getLinearCoordinate(int referenceIndex, long pos) {
		if (referenceIndex < 0 || referenceIndex >= offset.length) {
			return -1;
		}
		return offset[referenceIndex] + pos;
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getLinearCoordinate(java.lang.String, long)
	 */
	@Override
	public long getLinearCoordinate(String chr, long pos) {
		return getLinearCoordinate(dictionary.getSequenceIndex(chr), pos);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getStartLinearCoordinate(au.edu.wehi.idsv.BreakendSummary)
	 */
	@Override
	public long getStartLinearCoordinate(BreakendSummary bp) {
		return getLinearCoordinate(bp.referenceIndex, bp.start);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getStartLinearCoordinate(htsjdk.samtools.SAMRecord)
	 */
	@Override
	public long getStartLinearCoordinate(SAMRecord r) {
		return getLinearCoordinate(r.getReferenceIndex(), r.getAlignmentStart());
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getEndLinearCoordinate(au.edu.wehi.idsv.BreakendSummary)
	 */
	@Override
	public long getEndLinearCoordinate(BreakendSummary bp) {
		return getLinearCoordinate(bp.referenceIndex, bp.end);
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getEndLinearCoordinate(htsjdk.samtools.SAMRecord)
	 */
	@Override
	public long getEndLinearCoordinate(SAMRecord r) {
		return getLinearCoordinate(r.getReferenceIndex(), r.getAlignmentEnd());
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getReferenceIndex(long)
	 */
	@Override
	public int getReferenceIndex(long linearCoordinate) {
		return referenceIndices.floorEntry(linearCoordinate - 1).getValue();
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#getReferencePosition(long)
	 */
	@Override
	public int getReferencePosition(long linearCoordinate) {
		int referenceIndex = getReferenceIndex(linearCoordinate);
		if (referenceIndex < 0) return -1;
		return (int)(linearCoordinate - getLinearCoordinate(referenceIndex, 0));
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#encodedToString(long)
	 */
	@Override
	public String encodedToString(long linearCoordinate) {
		return String.format("%s:%d",
				dictionary.getSequence(getReferenceIndex(linearCoordinate)).getSequenceName(),
				getReferencePosition(linearCoordinate));
	}
	/* (non-Javadoc)
	 * @see au.edu.wehi.idsv.LinearGenomicCoordinate#encodedIntervalToString(long, long)
	 */
	@Override
	public String encodedIntervalToString(long startLinearCoordinate, long endLinearCoordinate) {
		if (startLinearCoordinate > endLinearCoordinate) throw new IllegalArgumentException("end cannot be before start");
		if (getReferenceIndex(startLinearCoordinate) == getReferenceIndex(endLinearCoordinate)) {
			return String.format("%s:%d-%d",
					dictionary.getSequence(getReferenceIndex(startLinearCoordinate)).getSequenceName(),
					getReferencePosition(startLinearCoordinate),
					getReferencePosition(endLinearCoordinate)); 
		} else {
			return encodedToString(startLinearCoordinate) + "-" + encodedToString(endLinearCoordinate);
		}
	}
}
