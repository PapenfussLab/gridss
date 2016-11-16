package au.edu.wehi.idsv;

import java.util.Iterator;
import java.util.List;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.configuration.RealignmentConfiguration;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * Incorporates realignments of soft clips into evidence 
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedSoftClipEvidenceIterator extends AbstractIterator<SoftClipEvidence> implements CloseableIterator<SoftClipEvidence>, TrackedBuffer {
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
		SoftClipEvidence evidence = breakendEvidence.next();
		if (evidence instanceof SpannedIndelEvidence) return evidence;
		RealignmentConfiguration rp = evidence.getEvidenceSource().getContext().getRealignmentParameters();
		SAMRecord realigned = factory.findFirstAssociatedSAMRecord(evidence, rp.shouldRealignBreakend(evidence));
		if (realigned != null) {
			return SoftClipEvidence.create(evidence, realigned);
		} else {
			return evidence;
		}
	}
	@Override
	public void close() {
		CloserUtil.close(breakendEvidence);
		CloserUtil.close(realigned);
	}
	@Override
	public void setTrackedBufferContext(String context) {
		factory.setTrackedBufferContext(context);
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return factory.currentTrackedBufferSizes();
	}
}
