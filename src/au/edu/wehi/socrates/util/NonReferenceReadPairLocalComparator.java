package au.edu.wehi.socrates.util;

import java.util.Comparator;

import au.edu.wehi.socrates.NonReferenceReadPair;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;

/**
 * Comparator for sorting NonReferenceReadPairs according to the local read
 * @author Daniel Cameron
 */
public class NonReferenceReadPairLocalComparator implements Comparator<NonReferenceReadPair> {
	private final Comparator<SAMRecord> comparator;
	public NonReferenceReadPairLocalComparator(Comparator<SAMRecord> comparator) {
		this.comparator = comparator;
	}
	@Override
	public int compare(NonReferenceReadPair arg0, NonReferenceReadPair arg1) {
		return comparator.compare(arg0.getLocalledMappedRead(), arg1.getLocalledMappedRead());
	}
}
