package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.Queue;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SplitIndel;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * Iterators over soft clip evidence in the same order as the input iterator 
 * @author Daniel Cameron
 *
 */
public class SoftClipEvidenceIterator implements CloseableIterator<SoftClipEvidence> {
	private final SAMEvidenceSource source;
	private final Iterator<SAMRecord> it;
	private final Queue<SoftClipEvidence> buffer = new ArrayDeque<SoftClipEvidence>(2);
	public SoftClipEvidenceIterator(SAMEvidenceSource source, Iterator<SAMRecord> it) {
		this.source = source;
		this.it = it;
	}
	private void fillBuffer() {
		while (it.hasNext() && buffer.isEmpty()) {
			SAMRecord record = it.next();
			if (!record.getReadUnmappedFlag()) {
				if (SAMRecordUtil.getStartSoftClipLength(record) > 0) {
					SoftClipEvidence sce = SoftClipEvidence.create(source, BreakendDirection.Backward, record);
					if (sce.meetsEvidenceCritera()) {
						buffer.add(sce);
					}
				}
				if (SAMRecordUtil.getEndSoftClipLength(record) > 0) {
					SoftClipEvidence sce = SoftClipEvidence.create(source, BreakendDirection.Forward, record);
					if (sce.meetsEvidenceCritera()) {
						buffer.add(sce);
					}
				} 
				int i = 0;
				for (SplitIndel si : SplitIndel.getIndelsAsSplitReads(record)) {
					SpannedIndelEvidence sie = SpannedIndelEvidence.create(source, record, si, i);
					buffer.add(sie);
					buffer.add(sie.asRemote());
					i++;
				}
			}
		}
	}
	@Override
	public boolean hasNext() {
		fillBuffer();
		return !buffer.isEmpty();
	}
	@Override
	public SoftClipEvidence next() {
		fillBuffer();
		return buffer.poll();
	}
	@Override
	public void close() {
		CloserUtil.close(it);
	}
}
