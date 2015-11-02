package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;

import au.edu.wehi.idsv.visualisation.TrackedBuffer;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

/**
 * Iterators through soft clip split reads in order of the remapped soft clipped position 
 * 
 * @author Daniel Cameron
 *
 */
public class RealignedRemoteSoftClipEvidenceIterator extends AbstractIterator<RealignedRemoteSoftClipEvidence> implements CloseableIterator<RealignedRemoteSoftClipEvidence>, TrackedBuffer {
	private final SAMEvidenceSource source;
	private final Iterator<SAMRecord> realigned;
	private final SequentialSoftClipRealignedRemoteBreakpointFactory ffactory;
	private final SequentialSoftClipRealignedRemoteBreakpointFactory bfactory;
	private final List<Object> toClose;
	private final Queue<RealignedRemoteSoftClipEvidence> buffer = new ArrayDeque<RealignedRemoteSoftClipEvidence>();
	public RealignedRemoteSoftClipEvidenceIterator(
			SAMEvidenceSource source,
			Iterator<SAMRecord> realigned,
			Iterator<SAMRecord> forwardSoftClipsOrderByRealignedCoordinate,
			Iterator<SAMRecord> backwardSoftClipsOrderByRealignedCoordinate) {
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
			String evidenceId = BreakpointFastqEncoding.getEncodedID(realign.getReadName());
			if (evidenceId.charAt(0) == BreakendDirection.Forward.toChar()) {
				RealignedRemoteSoftClipEvidence fsce = createRemote(BreakendDirection.Forward, ffactory.findFirstAssociatedSAMRecord(realign), realign);
				if (fsce != null) {
					buffer.add(fsce);
				}
			} else if (evidenceId.charAt(0) == BreakendDirection.Backward.toChar()) {
				RealignedRemoteSoftClipEvidence bsce = createRemote(BreakendDirection.Backward, bfactory.findFirstAssociatedSAMRecord(realign), realign);
				if (bsce != null) {
					buffer.add(bsce);
				}
			} else {
				throw new RuntimeException(String.format("Malformed realigned soft clip read name %s", evidenceId));
			}
		}
		if (buffer.isEmpty()) return endOfData();
		return buffer.poll();
	}
	private RealignedRemoteSoftClipEvidence createRemote(BreakendDirection direction, SAMRecord sc, SAMRecord realign) {
		if (sc == null || realign == null) return null;
		SoftClipEvidence breakendEvidence = SoftClipEvidence.create(source, direction, sc);
		if (breakendEvidence.meetsEvidenceCritera()) {
			SoftClipEvidence breakpointEvidence = SoftClipEvidence.create(breakendEvidence, realign);
			if (breakpointEvidence instanceof RealignedSoftClipEvidence) {
				return ((RealignedSoftClipEvidence)breakpointEvidence).asRemote();
			}
		}
		return null;
	}
	@Override
	public void close() {
		CloserUtil.close(toClose);
	}
	@Override
	public void setTrackedBufferContext(String context) {
		ffactory.setTrackedBufferContext(context + ".rsc.f");
		bfactory.setTrackedBufferContext(context + ".rsc.b");
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		List<NamedTrackedBuffer> list = new ArrayList<TrackedBuffer.NamedTrackedBuffer>();
		list.addAll(ffactory.currentTrackedBufferSizes());
		list.addAll(bfactory.currentTrackedBufferSizes());
		return list;
	}
}
