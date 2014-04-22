package au.edu.wehi.socrates;

import java.util.Comparator;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordComparator;
import net.sf.samtools.SAMRecordCoordinateComparator;

/**
 * Comparator for ordering @see DirectedEvidence by start position
 *
 */
public class DirectedEvidenceStartCoordinateComparator implements Comparator<DirectedEvidence> {
	@Override
	public int compare(DirectedEvidence arg0, DirectedEvidence arg1) {
		if (arg0.getReferenceIndex() < arg1.getReferenceIndex()) return -1;
		if (arg0.getReferenceIndex() == arg1.getReferenceIndex()) return (int)Math.abs(arg1.getWindowStart() - arg0.getWindowStart());
		return 1;
	}
}
