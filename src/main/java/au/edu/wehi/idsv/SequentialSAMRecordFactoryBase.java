package au.edu.wehi.idsv;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

/**
 * Common base class for supplying matching SAMRecords from SAMRecord sequence in
 * genomic coordinate order on non-standard key
 * @author Daniel Cameron
 * @param <T> type of records to find matching SAMRecord for. Must be suitable for placement in a @see HashMap
 * @param <U> type of matching key
 *
 */
public abstract class SequentialSAMRecordFactoryBase<T> implements TrackedBuffer {
	private static final Log log = Log.getInstance(SequentialSAMRecordFactoryBase.class);
	private final PeekingIterator<SAMRecord> sequence;
	private final HashMultimap<String, SAMRecord> currentReads = HashMultimap.create();
	private int currentReferenceIndex = -1; 
	private long currentPosition = -1;
	protected SequentialSAMRecordFactoryBase(PeekingIterator<SAMRecord> sequence) {
		if (sequence == null) {
			sequence = Iterators.peekingIterator(Collections.<SAMRecord>emptyIterator());
		}
		this.sequence = sequence;
	}
	protected SAMRecord findFirstMatching(int referenceIndex, int position, String key) {
		Set<SAMRecord> records = findMatching(referenceIndex, position, key);
		if (records == null || records.isEmpty()) return null;
		return records.iterator().next();
	}
	/**
	 * 
	 * @param referenceIndex expected reference index (according to the non-standard sort order) of matching SAMRecord
	 * @param position expected position (according to the non-standard sort order) of matching SAMRecord
	 * @param key key of matching SAMRecord
	 * @return matching SAMRecord if it exists, null otherwise
	 */
	protected Set<SAMRecord> findMatching(int referenceIndex, int position, String key) {
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
	protected abstract String getMatchingKey(SAMRecord record);
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
			int lastReferenceIndex = currentReferenceIndex;
			long lastPosition = currentPosition;
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
					checkSorted(lastReferenceIndex, lastPosition, currentReferenceIndex, currentPosition, r.getReadName());
					continue;
				} else if (sourceReferenceIndex == currentReferenceIndex) {
					if (sourcePosition < currentPosition) {
						// skip reads that are before our position
						checkSorted(lastReferenceIndex, lastPosition, currentReferenceIndex, currentPosition, r.getReadName());
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
	private void checkSorted(int lastReferenceIndex, long lastPosition, int referenceIndex, long position, String evidenceID) {
		if (referenceIndex < lastReferenceIndex || (referenceIndex == lastReferenceIndex && position < lastPosition)) {
			String msg = String.format("Unable to find realignment record for %s. This is likely due to either " 
					+ "a) alignment not completed successfully or "
					+ "b) chosen aligner writing records out of order."
					+ "The aligner is required to write records in same order as the input fastq. "
					+ "Please raise a github enhancement request if use by an aligner with no capability to ensure this ordering is required. ", evidenceID);
			log.equals(msg);
			throw new RuntimeException(msg);
		}
	}
	private String trackedBufferName = trackedName();
	protected abstract String trackedName();
	@Override
	public void setTrackedBufferContext(String context) {
		this.trackedBufferName = context + "." + trackedName();
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return ImmutableList.of(
				new NamedTrackedBuffer(trackedBufferName, currentReads.size())
				);
	}
}
