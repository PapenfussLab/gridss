package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class IdsvMetricsCollector {
	private IdsvMetrics metrics = new IdsvMetrics();
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
    	metrics.MAX_READ_LENGTH = Math.max(metrics.MAX_READ_LENGTH, record.getReadLength());
    	if (!record.getReadUnmappedFlag()) {
    		metrics.MAX_READ_MAPPED_LENGTH = Math.max(metrics.MAX_READ_MAPPED_LENGTH, record.getAlignmentEnd() - record.getAlignmentStart() + 1);
    	}
    	if (record.getReadPairedFlag() && record.getProperPairFlag()) {
    		int fragmentSize = SAMRecordUtil.estimateFragmentSize(record);
    		metrics.MAX_PROPER_PAIR_FRAGMENT_LENGTH = Math.max(metrics.MAX_PROPER_PAIR_FRAGMENT_LENGTH, Math.abs(fragmentSize));
    	}
    }
    public void finish() { }
	public void addAllLevelsToFile(MetricsFile<IdsvMetrics, Integer> imcf) {
		imcf.addMetric(metrics);
	}
}
