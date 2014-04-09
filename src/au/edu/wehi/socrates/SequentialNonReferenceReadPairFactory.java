package au.edu.wehi.socrates;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMRecord;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Iterates over read pairs ordered according to the start position of the locally mapped read in the pair
 * <p>Discordantly mapped read pairs will be iterated over twice, once of each read in the pair.</p>
 * @author Daniel Cameron
 *
 */
public class SequentialNonReferenceReadPairFactory {
	/**
	 * <p>read iterator <b>must</b> be coordinate sorted<p>
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param reads
	 * @param mate
	 */
	public SequentialNonReferenceReadPairFactory(
		Iterator<SAMRecord> mates) {
		this.mates = Iterators.peekingIterator(mates);
	}
	private PeekingIterator<SAMRecord> mates;
	public NonReferenceReadPair createNonReferenceReadPair(SAMRecord record) {
		if (!record.getReadPairedFlag()) return null;
		if (record.getProperPairFlag()) return null;
		NonReferenceReadPair pair = null;
		SAMRecord mate = findMate(record);
		if (record != null && mate != null) {
			return new NonReferenceReadPair(record, mate); 
		}
		return null;
	}
	private SAMRecord findMate(SAMRecord record) {
		if (record == null) return null;
		if (record.getMateUnmappedFlag()) return null;
		load(record.getMateReferenceIndex(), record.getMateAlignmentStart());
		for (SAMRecord mate : currentReads) {
			if (mate.getReadName() == record.getReadName() && mate.getFirstOfPairFlag() != record.getFirstOfPairFlag()) {
				return mate;
			}
		}
		return null;
	}
	private int currentReferenceIndex = -1; 
	private int currentPosition = -1;
	private List<SAMRecord> currentReads = new ArrayList<SAMRecord>();
	private void load(int referenceIndex, int alignmentStart) {
		if (referenceIndex > currentReferenceIndex || (referenceIndex == currentReferenceIndex && currentPosition > alignmentStart)) {
			throw new RuntimeException(String.format("Sanity check failure: cannot rewind NonReferenceReadPairIterator from %d:%d to %d:%d", currentReferenceIndex, currentPosition, referenceIndex, alignmentStart));
		}
		if (referenceIndex != currentReferenceIndex || currentPosition != alignmentStart) {
			currentReads.clear();
			currentReferenceIndex = referenceIndex;
			currentPosition = alignmentStart;
			// skip reads that are
			while (mates.hasNext() && (
					mates.peek().getMateUnmappedFlag() ||
					mates.peek().getMateReferenceIndex() < currentReferenceIndex ||
					(mates.peek().getMateReferenceIndex() == currentReferenceIndex && mates.peek().getMateAlignmentStart() == currentPosition))) {
				mates.next();
			}
			while (mates.hasNext() &&
					!mates.peek().getMateUnmappedFlag() &&
					mates.peek().getMateReferenceIndex() == currentReferenceIndex &&
					mates.peek().getMateAlignmentStart() == currentPosition) {
				currentReads.add(mates.next());
			}
		}
	}
}
