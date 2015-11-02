package au.edu.wehi.idsv;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.Queue;

import com.google.common.collect.Iterables;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

/**
 * Converts a ReadEvidenceAssembler to an iterator
 * 
 * @author Daniel Cameron
 *
 */
public class ReadEvidenceAssemblyIterator implements Iterator<SAMRecordAssemblyEvidence> {
	private static final Log log = Log.getInstance(ReadEvidenceAssemblyIterator.class);
	private static final int LOG_INTERVAL = 100000;
	private final ProgressLogger progress = new ProgressLogger(log);
	private final ReadEvidenceAssembler assembler;
	private final Iterator<DirectedEvidence> underlying;
	private final Queue<SAMRecordAssemblyEvidence> outputBuffer = new ArrayDeque<SAMRecordAssemblyEvidence>();
	private boolean flushed = false;
	public ReadEvidenceAssemblyIterator(ReadEvidenceAssembler assembler, Iterator<DirectedEvidence> evidence) {
		this.assembler = assembler;
		this.underlying = evidence;
	}
	private void ensureBuffer() {
		while (outputBuffer.isEmpty() && underlying.hasNext()) {
			Iterable<SAMRecordAssemblyEvidence> assemblies = assembler.addEvidence(track(underlying.next()));
			if (assemblies != null) {
				Iterables.addAll(outputBuffer, assemblies);
			}
		}
		if (outputBuffer.isEmpty() && !underlying.hasNext() && !flushed) {
			Iterable<SAMRecordAssemblyEvidence> assemblies = assembler.endOfEvidence();
			if (assemblies != null) {
				Iterables.addAll(outputBuffer, assemblies);
			}
			flushed = true;
		}
	}
	private DirectedEvidence track(DirectedEvidence readEvidence) {
		if (readEvidence instanceof NonReferenceReadPair) {
			track(((NonReferenceReadPair)readEvidence).getLocalledMappedRead());
		} else if (readEvidence instanceof SoftClipEvidence) {
			track(((SoftClipEvidence)readEvidence).getSAMRecord());
		}
		return readEvidence;
	}
	private long readCount = 0;
	private void track(SAMRecord record) {
		progress.record(record);
		readCount++;
		if (readCount % LOG_INTERVAL == 0) {
			log.info(String.format("After %d reads, assembly at %s:%d %s",
					readCount,
					record.getContig(), record.getAlignmentStart(),
					assembler.getStateSummaryMetrics()));
		}
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !outputBuffer.isEmpty();
	}
	@Override
	public SAMRecordAssemblyEvidence next() {
		ensureBuffer();
		return outputBuffer.poll();
	}
}
