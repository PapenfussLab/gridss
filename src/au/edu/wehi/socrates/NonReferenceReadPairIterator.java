package au.edu.wehi.socrates;

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
public class NonReferenceReadPairIterator extends AbstractIterator<NonReferenceReadPair> {
	/**
	 * <p>oea and dp iterators <b>must</b> be coordinate sorted<p>
	 * <p>mate iterator <b>must</b> be mate-coordinate sorted. @see SAMRecordMateCoordinateComparator<p>
	 * @param oea
	 * @param dp
	 * @param mate
	 */
	public NonReferenceReadPairIterator(
		Iterator<SAMRecord> oea,
		Iterator<SAMRecord> dp,
		Iterator<SAMRecord> mate) {
		this.reads = mergingIterator(oea, dp);
		this.mates = mate;
	}
	private Iterator<SAMRecord> reads;
	private Iterator<SAMRecord> mates;
	@Override
	protected NonReferenceReadPair computeNext() {
		if (!reads.hasNext()) return endOfData();
		SAMRecord next = reads.next();
		SAMRecord mate = findMate(next);
		return new NonReferenceReadPair(next, mate);
	}
	private SAMRecord findMate(SAMRecord record) {
		flushBefore(record.getReferenceIndex(), record.getAlignmentStart());
		return null;
	}
}
