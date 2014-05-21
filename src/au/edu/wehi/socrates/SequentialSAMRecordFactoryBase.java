package au.edu.wehi.socrates;

import htsjdk.samtools.SAMRecord;

import java.util.Map;

import com.google.common.collect.Iterators;
import com.google.common.collect.Maps;
import com.google.common.collect.PeekingIterator;

/**
 * Common base class for supplying matching SAMRecords from SAMRecord sequence in
 * genomic coordinate order on non-standard key
 * @author Daniel Cameron
 * @param <T> type of records to find matching SAMRecord for. Must be suitable for placement in a @see HashMap
 * @param <U> type of matching key
 *
 */
public abstract class SequentialSAMRecordFactoryBase<T> {
	private final PeekingIterator<SAMRecord> sequence;
	private final Map<T, SAMRecord> currentReads = Maps.newHashMap();
	private int currentReferenceIndex = -1; 
	private long currentPosition = -1;
	protected SequentialSAMRecordFactoryBase(
			PeekingIterator<SAMRecord> sequence) {
		if (sequence == null) {
			sequence = Iterators.peekingIterator(Iterators.<SAMRecord>emptyIterator());
		}
		this.sequence = sequence;
	}
	/**
	 * 
	 * @param referenceIndex expected reference index (according to the non-standard sort order) of matching SAMRecord
	 * @param position expected position (according to the non-standard sort order) of matching SAMRecord
	 * @param key key of matching SAMRecord
	 * @return matching SAMRecord if it exists, null otherwise
	 */
	protected SAMRecord findMatching(int referenceIndex, int position, T key) {
		load(referenceIndex, position);
		if (currentReads.containsKey(key)) {
			return currentReads.get(key);
		}
		return null;
	}
	/**
	 * Matching key on the sequence SAMRecord
	 * @param record SAMRecord to extract key from
	 * @return unique matching key for the given SAMRecord
	 */
	protected abstract T getMatchingKey(SAMRecord record);
	/**
	 * Non-standard genomic SAM contig index
	 * @param record
	 * @return reference index according to the non-standard coordinate sort order
	 */
	protected abstract int getReferenceIndex(SAMRecord record);
	/**
	 * Non-standard genomic position on contig
	 * @param record
	 * @return contig position according to the non-standard coordinate sort order
	 */
	protected abstract int getPosition(SAMRecord record);
	private void load(int referenceIndex, int alignmentStart) {
		if (referenceIndex < currentReferenceIndex || (referenceIndex == currentReferenceIndex && currentPosition > alignmentStart)) {
			throw new IllegalStateException(String.format("Sanity check failure: cannot rewind source iterator from %d:%d to %d:%d", currentReferenceIndex, currentPosition, referenceIndex, alignmentStart));
		}
		if (referenceIndex != currentReferenceIndex || currentPosition != alignmentStart) {
			currentReads.clear();
			currentReferenceIndex = referenceIndex;
			currentPosition = alignmentStart;
			while (sequence.hasNext()) {
				SAMRecord r = sequence.peek();
				int sourceReferenceIndex = getReferenceIndex(r);
				long sourcePosition = getPosition(r);
				if (sourceReferenceIndex < currentReferenceIndex) {
					// skip reads that are on earlier chromosome
					sequence.next();
					continue;
				} else if (sourceReferenceIndex == currentReferenceIndex) {
					if (sourcePosition < currentPosition) {
						// skip reads that are before our position
						sequence.next();
						continue;
					} else if (sourcePosition == currentPosition) {
						// add reads at current position
						currentReads.put(getMatchingKey(r), sequence.next());
						continue;
					} else {
						// stop when we hit a read after our position
						break;
					}
				} else {
					// stop when we hit a read on a subsequent chromosome
					break;
				}
			}
		}
	}
}
