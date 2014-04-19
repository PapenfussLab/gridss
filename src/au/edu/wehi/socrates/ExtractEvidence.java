package au.edu.wehi.socrates;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import au.edu.wehi.socrates.util.SAMRecordMateCoordinateComparator;
import au.edu.wehi.socrates.util.SAMRecordSummary;

import net.sf.picard.analysis.CollectInsertSizeMetrics;
import net.sf.picard.analysis.InsertSizeMetrics;
import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.analysis.directed.InsertSizeMetricsCollector;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CollectionUtil;
import net.sf.samtools.util.SortingCollection;

/**
 * Extracts reads supporting structural variation into intermediate files.
 * By sorting DP and OEA read pairs by the coordinate of the mapped mate,
 * putative directed breakpoints can be assembled by downstream processes
 * in a single pass of the intermediate files.
 * @author Daniel Cameron
 * @see GenerateDirectedBreakpoints
 *
 */
public class ExtractEvidence extends CommandLineProgram {
    private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Extract reads supporting structural variation." + PROGRAM_VERSION;

    @Option(doc="Coordinate-sorted input BAM file.",
    		shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    /*
    @Option(doc="BAM file containing reads supporting any putative structural variation breakpoint.",
    		shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;
    @Option(doc="Read pairs not mapped concordantly, sorted by mate coordinate.",
    		shortName="MCO")
    public File MATE_COORDINATE_OUTPUT;
    */
    @Option(doc = "Extract supporting reads into a separate file for each chromosome.",
            optional = true,
            shortName = "BYCHR")
    public boolean PER_CHR = true;
    private Log log = Log.getInstance(ExtractEvidence.class);
    private void writeToFile(SortingCollection<SAMRecord> sc, SAMFileWriter writer) {
    	sc.doneAdding();
    	final CloseableIterator<SAMRecord> it = sc.iterator();
		while (it.hasNext()) {
			writer.addAlignment(it.next());
		}
    }
    @Override
	protected int doWork() {
    	try {
	    	IoUtil.assertFileIsWritable(FileNamingConvention.GetMetrics(INPUT));
	    	SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	    	// SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(CREATE_INDEX); // do we need an index at all? Can we get away with sequential access?
	    	final SAMFileReader reader = new SAMFileReader(INPUT);
	    	final SAMFileHeader header = reader.getFileHeader();
	    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
	    	final InsertSizeMetricsCollector metrics = RelevantMetrics.createCollector(header);
	    	
	    	final SAMFileWriter[] writers = new SAMFileWriter[dictionary.size() + 1];
	    	final SAMFileWriter[] matewriters = new SAMFileWriter[dictionary.size() + 1];
	    	final ArrayList<SortingCollection<SAMRecord>> matecollection = new ArrayList<SortingCollection<SAMRecord>>();
	    	final SAMFileHeader svHeader = header.clone();
	    	svHeader.setSortOrder(SortOrder.coordinate);
	    	final SAMFileHeader mateHeader = header.clone();
	    	mateHeader.setSortOrder(SortOrder.unsorted);
	    	final SAMFileWriterFactory factory = new SAMFileWriterFactory();
	    	if (PER_CHR) {
	    		for (int i = 0; i < dictionary.size(); i++) {
	    			writers[i] = factory.makeBAMWriter(svHeader, true, FileNamingConvention.GetSVBamForChr(INPUT, dictionary.getSequence(i).getSequenceName()));
	    			matewriters[i] = factory.makeBAMWriter(mateHeader, false, FileNamingConvention.GetMateBamForChr(INPUT, dictionary.getSequence(i).getSequenceName()));
	    			matecollection.add(SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(mateHeader), new SAMRecordMateCoordinateComparator(),
	    					// TODO: allocate buffers according to sequence lengths instead of equal space to every chr
	    					MAX_RECORDS_IN_RAM / dictionary.size(),
	    					TMP_DIR));
	    		}
				writers[dictionary.size()] = factory.makeBAMWriter(svHeader, true, FileNamingConvention.GetSVBamForChr(INPUT, "unmapped"));
	    	} else {
	    		// all writers map to the same one
	    		writers[0] = factory.makeBAMWriter(svHeader, true, FileNamingConvention.GetSVBam(INPUT));
				matewriters[0] = factory.makeBAMWriter(mateHeader, false, FileNamingConvention.GetMateBam(INPUT));
				matecollection.add(SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(matewriters[0].getFileHeader()), new SAMRecordMateCoordinateComparator(), MAX_RECORDS_IN_RAM, TMP_DIR));
				for (int i = 1; i < dictionary.size() + 1; i++) {
					writers[i] = writers[0];
					matewriters[i] = matewriters[0];
				}
	    	}
	    	// Traverse the input file
	    	final ProgressLogger progress = new ProgressLogger(log);
			final SAMRecordIterator iter = reader.iterator();
			iter.assertSorted(SortOrder.coordinate);
			while (iter.hasNext()) {
				SAMRecord record = iter.next();
				int offset = record.getReadUnmappedFlag() ? dictionary.size() : record.getReferenceIndex();
				boolean sc = SAMRecordSummary.isAlignmentSoftClipped(record);
				boolean badpair = SAMRecordSummary.isPartOfNonReferenceReadPair(record);
				if (sc || badpair) {
					writers[offset].addAlignment(record);
				}
				if (badpair) {
					matecollection.get(record.getMateReferenceIndex()).add(record);
				}
				metrics.acceptRecord(record, null);
				progress.record(record);
			}
			reader.close();
			metrics.finish();
			RelevantMetrics.save(metrics, this.<InsertSizeMetrics, Integer>getMetricsFile(), FileNamingConvention.GetMetrics(INPUT));
			for (int i = 0; i < dictionary.size(); i++) {
				if (writers[i] != null) {
					writers[i].close();
					writers[i] = null;
				}
			}
			if (PER_CHR) {
				for (int i = 0; i < dictionary.size(); i++) {
					// TODO: this can be multi-threaded
					matewriters[i].setProgressLogger(new ProgressLogger(log));
					writeToFile(matecollection.get(i), matewriters[i]);
					matewriters[i].close();
					matecollection.set(i, null);
				}
				
			} else {
				matewriters[0].setProgressLogger(new ProgressLogger(log));
				writeToFile(matecollection.get(0), matewriters[0]);
				matewriters[0].close();
			}
			matecollection.clear();
    	} catch (IOException e) {
    		throw new RuntimeException(e);
    	}
        return 0;
    }
	public static void main(String[] argv) {
        System.exit(new ExtractEvidence().instanceMain(argv));
    }
}
