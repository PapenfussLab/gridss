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
	public final SAMRecord rightRealigned;
	public final SAMRecord rightAnchored;
	public final SAMRecord leftRealigned;
	public SplitIndel(SAMRecord read, int leftAnchorOffset, int rightAnchorOffset) {
		Pair<SAMRecord, SAMRecord> left = SAMRecordUtil.splitAfter(read, leftAnchorOffset);
		this.leftAnchored = left.getLeft();
		this.rightRealigned = left.getRight();
		SAMRecordUtil.trim(this.rightRealigned, this.leftAnchored.getReadLength() - SAMRecordUtil.getEndSoftClipLength(this.leftAnchored), 0);
		
		Pair<SAMRecord, SAMRecord> right = SAMRecordUtil.splitAfter(read, rightAnchorOffset - 1);
		this.rightAnchored = right.getLeft();
		this.leftRealigned = right.getRight();
		SAMRecordUtil.trim(this.leftRealigned, 0, this.rightAnchored.getReadLength() - SAMRecordUtil.getStartSoftClipLength(this.rightAnchored));
	}
	/**
	 * Converts all CIGAR indels into split indels
	 * @return
	 */
	public static List<SplitIndel> getIndelsAsSplitReads(SAMRecord read) {
		List<SplitIndel> list = new ArrayList<SplitIndel>(CigarUtil.countIndels(read.getCigar()));
		int offset = 0;
		for (CigarElement ce : CigarUtil.decodeNegativeDeletion(read.getCigar().getCigarElements())) {
			if (ce.getOperator() == CigarOperator.DELETION || ce.getOperator() == CigarOperator.INSERTION) {
				list.add(new SplitIndel(read, offset, ce.getOperator().consumesReadBases() ? offset + ce.getLength() + 1 : 1));
			}
			if (ce.getOperator().consumesReadBases()) {
				offset += ce.getLength();
			}
		}
		return list;
	}
}
