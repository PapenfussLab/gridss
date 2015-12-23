package au.edu.wehi.idsv.sam;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Converts read alignment indels into split reads
 * @author cameron.d
 *
 */
public class SplitIndel {
	public final SAMRecord leftAnchored;
	public final SAMRecord leftRealigned;
	public final SAMRecord rightAnchored;
	public final SAMRecord rightRealigned;
	public SplitIndel(SAMRecord read, int anchorOffset) {
		Pair<SAMRecord, SAMRecord> left = SAMRecordUtil.splitAfter(read, anchorOffset);
		this.leftAnchored = left.getLeft();
		this.leftRealigned = left.getRight();
		SAMRecordUtil.trim(this.leftRealigned, this.leftAnchored.getReadLength() - SAMRecordUtil.getEndSoftClipLength(this.leftAnchored), 0);
		Pair<SAMRecord, SAMRecord> right = SAMRecordUtil.splitAfter(read, anchorOffset);
		this.rightRealigned = right.getLeft();
		this.rightAnchored = right.getRight();
		SAMRecordUtil.trim(this.rightRealigned, 0, this.rightAnchored.getReadLength() - SAMRecordUtil.getStartSoftClipLength(this.rightAnchored));
		sanityCheck();
	}
	private void sanityCheck() {
		assert(leftAnchored.getAlignmentStart() == rightRealigned.getAlignmentStart());
		assert(rightAnchored.getAlignmentStart() == leftRealigned.getAlignmentStart());
		assert(leftAnchored.getReadLength() == rightAnchored.getReadLength());
	}
	/**
	 * Converts all CIGAR indels into split indels
	 * @return
	 */
	public static List<SplitIndel> getIndelsAsSplitReads(SAMRecord read) {
		List<SplitIndel> list = new ArrayList<SplitIndel>(CigarUtil.countIndels(read.getCigar()));
		int offset = 0;
		boolean inIndel = false;
		for (CigarElement ce : CigarUtil.decodeNegativeDeletion(read.getCigar().getCigarElements())) {
			if (ce.getOperator() == CigarOperator.DELETION || ce.getOperator() == CigarOperator.INSERTION) {
				if (!inIndel) {
					list.add(new SplitIndel(read, offset - 1));
				}
				inIndel = true;
			} else {
				inIndel = false;
			}
			if (ce.getOperator().consumesReadBases()) {
				offset += ce.getLength();
			}
		}
		return list;
	}
}
