package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;
import java.util.Queue;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Queues;

/**
 * Iterators over soft clip evidence
 * @author Daniel Cameron
 *
 */
public class SoftClipEvidenceIterator extends AbstractIterator<SoftClipEvidence> implements CloseableIterator<SoftClipEvidence> {
	private final ProcessingContext processContext;
	private final SAMEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final Queue<SoftClipEvidence> buffer = Queues.newArrayDeque();
	public SoftClipEvidenceIterator(ProcessingContext processContext, SAMEvidenceSource source, Iterator<SAMRecord> it) {
		this.processContext = processContext;
		this.source = source;
		this.it = it;
	}
	private void fillBuffer() {
		while (it.hasNext() && buffer.isEmpty()) {
			SAMRecord record = it.next();
			if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
				SoftClipEvidence sce = SoftClipEvidence.create(processContext, source, BreakendDirection.Backward, record);
				if (processContext.getSoftClipParameters().meetsEvidenceCritera(sce)) {
					buffer.add(sce);
				}
			}
			if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
				SoftClipEvidence sce = SoftClipEvidence.create(processContext, source, BreakendDirection.Forward, record);
				if (processContext.getSoftClipParameters().meetsEvidenceCritera(sce)) {
					buffer.add(sce);
				}
			}
		}
	}
	@Override
	protected SoftClipEvidence computeNext() {
		if (buffer.isEmpty()) fillBuffer();
		if (buffer.isEmpty()) return endOfData();
		return buffer.poll();
	}
	@Override
	public void close() {
		CloserUtil.close(it);
	}
}
