package au.edu.wehi.idsv.sam;

import java.util.Iterator;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import htsjdk.samtools.SAMRecord;

/**
 * Ensures that every aligned record has an associated NM tag 
 * 
 * @author Daniel Cameron
 *
 */
public class NmTagIterator implements Iterator<SAMRecord> {
	private final Iterator<SAMRecord> it;
	private final ReferenceLookup reference;
	public NmTagIterator(Iterator<SAMRecord> it, ReferenceLookup reference) {
		this.it = it;
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
