package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

/**
 * Iterators through soft clip split reads in order of the remapped soft clipped position 
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedRemoteSoftClipEvidenceIterator extends AbstractIterator<RealignedRemoteSoftClipEvidence> implements CloseableIterator<RealignedRemoteSoftClipEvidence> {
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final Iterator<SAMRecord> realigned;
	private final SequentialSoftClipRealignedRemoteBreakpointFactory ffactory;
	private final SequentialSoftClipRealignedRemoteBreakpointFactory bfactory;
	private final List<Object> toClose;
	private final Queue<RealignedRemoteSoftClipEvidence> buffer = new ArrayDeque<RealignedRemoteSoftClipEvidence>();
	public RealignedRemoteSoftClipEvidenceIterator(
			ProcessingContext processContext,
			SAMEvidenceSource source,
			Iterator<SAMRecord> realigned,
			Iterator<SAMRecord> forwardSoftClipsOrderByRealignedCoordinate,
			Iterator<SAMRecord> backwardSoftClipsOrderByRealignedCoordinate) {
		this.processContext = processContext;
		this.source = source;
		this.toClose = Lists.<Object>newArrayList(realigned, forwardSoftClipsOrderByRealignedCoordinate, backwardSoftClipsOrderByRealignedCoordinate);
		this.realigned = realigned;
		this.ffactory = new SequentialSoftClipRealignedRemoteBreakpointFactory(Iterators.peekingIterator(forwardSoftClipsOrderByRealignedCoordinate), BreakendDirection.Forward);
		this.bfactory = new SequentialSoftClipRealignedRemoteBreakpointFactory(Iterators.peekingIterator(backwardSoftClipsOrderByRealignedCoordinate), BreakendDirection.Backward);
	}
	@Override
	protected RealignedRemoteSoftClipEvidence computeNext() {
		while (realigned.hasNext() && buffer.isEmpty()) {
			SAMRecord realign = realigned.next();
			SAMRecord fsc = ffactory.findAssociatedSAMRecord(realign);
			if (fsc != null) {
				buffer.add(new RealignedRemoteSoftClipEvidence(processContext, source, BreakendDirection.Forward, fsc, realign));
			}
			SAMRecord bsc = bfactory.findAssociatedSAMRecord(realign);
			if (bsc != null) {
				buffer.add(new RealignedRemoteSoftClipEvidence(processContext, source, BreakendDirection.Backward, bsc, realign));
			}
		}
		if (buffer.isEmpty()) return endOfData();
		return buffer.poll();
	}
	@Override
	public void close() {
		CloserUtil.close(toClose);
	}
}
