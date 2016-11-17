package au.edu.wehi.idsv.metrics;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.google.common.collect.Iterables;

import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import gridss.analysis.CigarDetailMetrics;
import gridss.analysis.IdsvMetrics;
import gridss.analysis.MapqMetrics;
import gridss.analysis.MapqMetricsCollector;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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
 * @author Daniel Cameron
 *
 */
public class IdsvSamFileMetricsCollector {
	private IdsvMetrics idsv = new IdsvMetrics();
	private HashMap<CigarOperator, List<CigarDetailMetrics>> cigar;
	private InsertSizeMetricsCollector is;
	private MapqMetricsCollector mmc;
	public IdsvSamFileMetricsCollector(SAMFileHeader header) {
		this.is = createInsertSizeMetricsCollector(header);
		this.mmc = createMapqMetricsCollector(header);
		this.cigar = new HashMap<CigarOperator, List<CigarDetailMetrics>>();
		for (CigarOperator op : CigarOperator.values()) {
			cigar.put(op, new ArrayList<CigarDetailMetrics>());
		}
	}
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
    	is.acceptRecord(record, refSeq);
    	idsvAcceptRecord(record, refSeq);
    	cigarAcceptRecord(record, refSeq);
    	mapqAcceptRecord(record, refSeq);
    }
    private void mapqAcceptRecord(SAMRecord record, ReferenceSequence refSeq) {
		mmc.acceptRecord(record, refSeq);
	}
	private void cigarAcceptRecord(SAMRecord record, ReferenceSequence refSeq) {
    	if (record == null || record.getCigar() == null) return;
    	List<CigarElement> list = record.getCigar().getCigarElements();
    	if (list == null || list.size() == 0) return;
    	for (CigarElement ce : list) {
    		acceptCigarElement(ce);
    	}
    	for (CigarOperator op : CigarOperator.values()) {
    		switch (op) {
    			case S:
    				if (CigarUtil.getStartClipLength(list) == 0) {
    					acceptCigarElement(new CigarElement(0, CigarOperator.S));
    				}
    				if (CigarUtil.getEndClipLength(list) == 0) {
    					acceptCigarElement(new CigarElement(0, CigarOperator.S));
    				}
    				break;
    			case H:
    				if (list.get(0).getOperator() != CigarOperator.H) {
    					acceptCigarElement(new CigarElement(0, CigarOperator.H));
    				}
    				if (list.get(list.size() - 1).getOperator() != CigarOperator.H) {
    					acceptCigarElement(new CigarElement(0, CigarOperator.H));
    				}
    				break;
    			default:
    				if (!Iterables.any(list, ce -> ce.getOperator() == op)) {
    					acceptCigarElement(new CigarElement(0, op));
    				}
    				break;
    		}
    	}
	}
    private void acceptCigarElement(CigarElement ce) {
    	List<CigarDetailMetrics> list = cigar.get(ce.getOperator());
    	int length = ce.getLength();
    	while (list.size() <= length) {
    		CigarDetailMetrics cdm = new CigarDetailMetrics();
    		cdm.LENGTH = list.size();
    		cdm.OPERATOR = (char)CigarOperator.enumToCharacter(ce.getOperator());
    		cdm.COUNT = 0;
    		list.add(cdm);
    	}
    	list.get(ce.getLength()).COUNT++;
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
		MetricsFile<CigarDetailMetrics, Integer> scMetricsFile = processContext.<CigarDetailMetrics, Integer>createMetricsFile();
		MetricsFile<MapqMetrics, Integer> mapqMetricsFile = processContext.<MapqMetrics, Integer>createMetricsFile();
		
		finish(isMetricsFile, idsvMetricsFile, mapqMetricsFile, scMetricsFile);
		
		isMetricsFile.write(processContext.getFileSystemContext().getInsertSizeMetrics(source));
		idsvMetricsFile.write(processContext.getFileSystemContext().getIdsvMetrics(source));
		scMetricsFile.write(processContext.getFileSystemContext().getCigarMetrics(source));
		mapqMetricsFile.write(processContext.getFileSystemContext().getMapqMetrics(source));
	}
    public void finish(MetricsFile<InsertSizeMetrics, Integer> isMetricsFile, MetricsFile<IdsvMetrics, Integer> idsvMetricsFile, MetricsFile<MapqMetrics, Integer> mapqMetricsFile, MetricsFile<CigarDetailMetrics, Integer> scMetricsFile) {
    	addInsertSizeMetrics(isMetricsFile);
		addIdsvMetrics(idsvMetricsFile);
		addCigarMetrics(scMetricsFile);
		addMapqMetrics(mapqMetricsFile);
    }
    private void addMapqMetrics(MetricsFile<MapqMetrics, Integer> mapqMetricsFile) {
    	mmc.finish();
    	mmc.addAllLevelsToFile(mapqMetricsFile);
	}
	private void addInsertSizeMetrics(MetricsFile<InsertSizeMetrics, Integer> metricsFile) {
    	is.finish();
    	is.addAllLevelsToFile(metricsFile);
	}
	private void addIdsvMetrics(MetricsFile<IdsvMetrics, Integer> metricsFile) {
		metricsFile.addMetric(idsv);
	}
	private void addCigarMetrics(MetricsFile<CigarDetailMetrics, Integer> metricsFile) {
		cigar.values().stream().flatMap(c -> c.stream()).forEach(metric -> {
			metricsFile.addMetric(metric);
		});
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
				new CollectInsertSizeMetrics().HISTOGRAM_WIDTH,
				new CollectInsertSizeMetrics().DEVIATIONS,
				true);
	}
	private static MapqMetricsCollector createMapqMetricsCollector(SAMFileHeader header) {
		return new MapqMetricsCollector(CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS), null);
	}
}
