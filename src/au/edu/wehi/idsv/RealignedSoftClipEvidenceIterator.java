package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

/**
 * Incorporates realignments of soft clips into evidence 
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedSoftClipEvidenceIterator extends AbstractIterator<SoftClipEvidence> implements CloseableIterator<SoftClipEvidence> {
	private final SequentialRealignedBreakpointFactory factory;
	private final Iterator<SoftClipEvidence> breakendEvidence;
	private final Iterator<SAMRecord> realigned;
	/**
	 * Create new iterator
	 * @param breakendEvidence soft clip breakend evidence sorted by breakend genomic position
	 * @param realignedSoftClips realigned soft clips in breakend evidence order
	 */
	public RealignedSoftClipEvidenceIterator(
			Iterator<SoftClipEvidence> breakendEvidence,
			Iterator<SAMRecord> realignedSoftClips) {
		this.breakendEvidence = breakendEvidence;
		this.realigned = realignedSoftClips;
		this.factory = new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(realignedSoftClips));
	}
	@Override
	protected SoftClipEvidence computeNext() {
		if (!breakendEvidence.hasNext()) return endOfData();
		SoftClipEvidence sce = breakendEvidence.next();
		SAMRecord realigned = factory.findAssociatedSAMRecord(sce);
		if (realigned != null) {
			return SoftClipEvidence.create(sce, realigned);
		} else {
			return sce;
		}
	}
	@Override
	public void close() {
		CloserUtil.close(breakendEvidence);
		CloserUtil.close(realigned);
	}
}
