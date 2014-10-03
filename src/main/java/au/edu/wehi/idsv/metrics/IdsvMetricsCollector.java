package au.edu.wehi.idsv.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class IdsvMetricsCollector {
	private IdsvMetrics metrics = new IdsvMetrics();
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
    	metrics.MAX_READ_LENGTH = Math.max(metrics.MAX_READ_LENGTH, record.getReadLength());
    	if (record.getReadPairedFlag() && record.getProperPairFlag()) {
    		// getInferredInsertSize is 5' to 5' and only work for FR
    		// ideally we want the unclipped start/end of the mate but we don't have access to that
    		int fragmentSize;
    		if (record.getUnclippedStart() < record.getMateAlignmentStart()) {
    			// will overestimate when mate is soft clipped at start
    			fragmentSize = record.getMateAlignmentStart() + record.getReadLength() - record.getUnclippedStart();
    		} else {
    			// will underestimate when mate is soft clipped at start
    			fragmentSize = record.getUnclippedEnd() - record.getMateAlignmentStart();
    		}
    		metrics.MAX_PROPER_PAIR_FRAGMENT_LENGTH = Math.max(metrics.MAX_PROPER_PAIR_FRAGMENT_LENGTH, Math.abs(fragmentSize));
    	}
    }
    public void finish() { }
	public void addAllLevelsToFile(MetricsFile<IdsvMetrics, Integer> imcf) {
		imcf.addMetric(metrics);
	}
}
