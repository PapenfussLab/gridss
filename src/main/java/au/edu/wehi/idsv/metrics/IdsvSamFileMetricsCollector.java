package au.edu.wehi.idsv.metrics;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import picard.analysis.CollectInsertSizeMetrics;
import picard.analysis.InsertSizeMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.directed.InsertSizeMetricsCollector;

/**
 * Collects metrics required by gridss
 * 
 * @author cameron.d
 *
 */
public class IdsvSamFileMetricsCollector {
	private IdsvMetrics idsv = new IdsvMetrics();
	private List<SoftClipDetailMetrics> sc = new ArrayList<SoftClipDetailMetrics>();
	private InsertSizeMetricsCollector is;
	public IdsvSamFileMetricsCollector(SAMFileHeader header) {
		this.is = createInsertSizeMetricsCollector(header);
	}
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
    	is.acceptRecord(record, refSeq);
    	idsvAcceptRecord(record, refSeq);
    	scAcceptRecord(record, refSeq);
    }
    private void scAcceptRecord(SAMRecord record, ReferenceSequence refSeq) {
    	while (sc.size() < record.getReadLength()) {
    		SoftClipDetailMetrics scNext = new SoftClipDetailMetrics();
    		scNext.LENGTH = sc.size();
    		sc.add(scNext);
    	}
    	if (!record.getReadUnmappedFlag()) {
	    	sc.get(SAMRecordUtil.getEndSoftClipLength(record)).READCOUNT++;
	    	sc.get(SAMRecordUtil.getStartSoftClipLength(record)).READCOUNT++;
    	}
	}
	private void idsvAcceptRecord(SAMRecord record, ReferenceSequence refSeq) {
		idsv.MAX_READ_LENGTH = Math.max(idsv.MAX_READ_LENGTH, record.getReadLength());
    	if (!record.getReadUnmappedFlag()) {
    		idsv.MAX_READ_MAPPED_LENGTH = Math.max(idsv.MAX_READ_MAPPED_LENGTH, record.getAlignmentEnd() - record.getAlignmentStart() + 1);
    	}
    	if (record.getReadPairedFlag()) {
    		if (record.getProperPairFlag()) {
	    		int fragmentSize = SAMRecordUtil.estimateFragmentSize(record, PairOrientation.FR);
	    		fragmentSize = Math.abs(fragmentSize);
	    		if (idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH == null) {
	    			idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH = fragmentSize;
	    		} else {
	    			idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH = Math.max(idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH, Math.abs(fragmentSize));
	    		}
	    		if (idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH == null) {
	    			idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH = fragmentSize;
	    		} else {
	    			idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH = Math.min(idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH, Math.abs(fragmentSize));
	    		}
    		}
    		if (record.getFirstOfPairFlag()) {
    			idsv.READ_PAIRS++;
    			if (record.getReadUnmappedFlag() && record.getMateUnmappedFlag()) {
    				idsv.READ_PAIRS_ZERO_MAPPED++;
    			} else if (!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
    				idsv.READ_PAIRS_BOTH_MAPPED++;
    			} else {
    				idsv.READ_PAIRS_ONE_MAPPED++;
    			}
    		}
    	}
    	idsv.READS++;
    	if (!record.getReadUnmappedFlag()) {
    		idsv.MAPPED_READS++;
    	}
	}
	public void finish(ProcessingContext processContext, File source) {		
		MetricsFile<InsertSizeMetrics, Integer> isMetricsFile = processContext.<InsertSizeMetrics, Integer>createMetricsFile();
		MetricsFile<IdsvMetrics, Integer> idsvMetricsFile = processContext.<IdsvMetrics, Integer>createMetricsFile();
		MetricsFile<SoftClipDetailMetrics, Integer> scMetricsFile = processContext.<SoftClipDetailMetrics, Integer>createMetricsFile();
		
		finish(isMetricsFile, idsvMetricsFile, scMetricsFile);
		
		isMetricsFile.write(processContext.getFileSystemContext().getInsertSizeMetrics(source));
		idsvMetricsFile.write(processContext.getFileSystemContext().getIdsvMetrics(source));
		scMetricsFile.write(processContext.getFileSystemContext().getSoftClipMetrics(source));
	}
    public void finish(MetricsFile<InsertSizeMetrics, Integer> isMetricsFile, MetricsFile<IdsvMetrics, Integer> idsvMetricsFile, MetricsFile<SoftClipDetailMetrics, Integer> scMetricsFile) {
    	addInsertSizeMetrics(isMetricsFile);
		addIdsvMetrics(idsvMetricsFile);
		addSoftClipMetrics(scMetricsFile);
    }
    private void addInsertSizeMetrics(MetricsFile<InsertSizeMetrics, Integer> metricsFile) {
    	is.finish();
    	is.addAllLevelsToFile(metricsFile);
	}
	private void addIdsvMetrics(MetricsFile<IdsvMetrics, Integer> metricsFile) {
		metricsFile.addMetric(idsv);
	}
	private void addSoftClipMetrics(MetricsFile<SoftClipDetailMetrics, Integer> metricsFile) {
		for (SoftClipDetailMetrics metrics : sc)  {
			metricsFile.addMetric(metrics);
		}
	}
	private static InsertSizeMetricsCollector createInsertSizeMetricsCollector(SAMFileHeader header) {
		//List<SAMReadGroupRecord> rg = ImmutableList.<SAMReadGroupRecord>of();
		//if (header != null) {
		//	rg = header.getReadGroups();
		//}
		return new InsertSizeMetricsCollector(
    			CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS), null, //, MetricAccumulationLevel.SAMPLE), rg,
				// match CollectInsertSizeMetrics defaults
				new CollectInsertSizeMetrics().MINIMUM_PCT,
				new CollectInsertSizeMetrics().Histogram_WIDTH,
				new CollectInsertSizeMetrics().DEVIATIONS);
	}
}
