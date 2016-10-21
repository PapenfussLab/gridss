package au.edu.wehi.idsv.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;

import java.util.Iterator;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Ensures that every aligned record has an associated NM tag 
 * 
 * @author Daniel Cameron
 *
 */
public class NmTagIterator implements Iterator<SAMRecord> {
	private final PeekingIterator<SAMRecord> it;
	private final ReferenceSequenceFile reference;
	public NmTagIterator(Iterator<SAMRecord> it, ReferenceSequenceFile reference) {
		this.it = Iterators.peekingIterator(it);
		this.reference = reference;
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}
	@Override
	public SAMRecord next() {
		SAMRecord r = it.next();
		SAMRecordUtil.ensureNmTag(reference, r);
		return r;
	}

}
