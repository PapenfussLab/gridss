package au.edu.wehi.socrates;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import javax.sound.midi.Sequence;

import au.edu.wehi.socrates.FileNamingConvention.IntermediateBamType;
import au.edu.wehi.socrates.util.SAMFileInfo;
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
import net.sf.picard.filter.AggregateFilter;
import net.sf.picard.filter.DuplicateReadFilter;
import net.sf.picard.filter.FailsVendorReadQualityFilter;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMException;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordQueryNameComparator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CollectionUtil;
import net.sf.samtools.util.IOUtil;
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

    @Option(doc="The BAM file to process.",
    		shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    /*// Output file names determine by @see FileNamingConvention
    @Option(doc = "Coordinate sorted output file for soft clipped reads",
            optional = false,
            shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File SC_OUTPUT;
    @Option(doc = "Coordinate sorted output file for discordantly paired reads",
            optional = false,
            shortName = "DPO")
    public File DP_OUTPUT = null;
    @Option(doc = "Coordinate sorted output file for reads pairs with only a single read mapped (open ended anchor)",
            optional = false,
            shortName = "OEAO")
    public File OEA_OUTPUT = null;
	@Option(doc = "DP and OEA reads sorted by coordinate of mapped mate read.",
            optional = true,
            shortName = "MCO")
    public File MATE_COORDINATE_OUTPUT = null;
    */
    @Option(doc = "Extract supporting reads into a separate file for each chromosome.",
            optional = true,
            shortName = "BYCHR")
    public boolean PER_CHR = true;
    @Option(doc = "Extract supporting reads into a separate file for each read group.",
            optional = true,
            shortName = "BYRG")
    public boolean PER_READ_GROUP = false;
    private Log log = Log.getInstance(ExtractEvidence.class);
    private SAMFileWriter getCoordinateSortedWriter(SAMFileHeader header, File file) {
    	SAMFileHeader outputHeader = header.clone();
    	outputHeader.setSortOrder(SortOrder.coordinate);
	    final SAMFileWriter writer = new SAMFileWriterFactory()
	    	.makeBAMWriter(outputHeader, true, file);
	    return writer;
    }
    private SAMFileWriter getMateCoordinateSortedWriter(SAMFileHeader header, File file) {
    	SAMFileHeader mateHeader = header.clone();
		mateHeader.setSortOrder(SortOrder.unsorted);
	    final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(mateHeader, false, file);
	    return writer;
    }
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
	    	IoUtil.assertFileIsWritable(FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.SC));
	    	IoUtil.assertFileIsWritable(FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.OEA));
	    	IoUtil.assertFileIsWritable(FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.DP));
	    	IoUtil.assertFileIsWritable(FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.MATE));
	    	if (PER_READ_GROUP) {
	    		throw new RuntimeException("Not Yet Implemented");
	    	}
	    	SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	    	// SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(CREATE_INDEX); // do we need an index at all? Can we get away with sequential access?
	    	final SAMFileReader reader = new SAMFileReader(INPUT);
	    	final SAMFileHeader header = reader.getFileHeader();
	    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
	    	final InsertSizeMetricsCollector metrics = new InsertSizeMetricsCollector(
	    			CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS, MetricAccumulationLevel.READ_GROUP),
					header.getReadGroups(),
					// match CollectInsertSizeMetrics defaults
					 new CollectInsertSizeMetrics().MINIMUM_PCT,
					 new CollectInsertSizeMetrics().HISTOGRAM_WIDTH,
					 new CollectInsertSizeMetrics().DEVIATIONS);
	    	
	    	SAMFileWriter[] scwriters = new SAMFileWriter[dictionary.size() + 1];
	    	SAMFileWriter[] oeawriters = new SAMFileWriter[dictionary.size() + 1];
	    	SAMFileWriter[] dpwriters = new SAMFileWriter[dictionary.size() + 1];
	    	SAMFileWriter[] matewriters = new SAMFileWriter[dictionary.size() + 1];
	    	ArrayList<SortingCollection<SAMRecord>> matecollection = new ArrayList<SortingCollection<SAMRecord>>();
	    	if (PER_CHR) {
	    		for (int i = 0; i < dictionary.size(); i++) {
	    			scwriters[i] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBamForChr(INPUT, IntermediateBamType.SC, dictionary.getSequence(i).getSequenceName()));
	    			oeawriters[i] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBamForChr(INPUT, IntermediateBamType.OEA, dictionary.getSequence(i).getSequenceName()));
	    			dpwriters[i] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBamForChr(INPUT, IntermediateBamType.DP, dictionary.getSequence(i).getSequenceName()));
	    			matewriters[i] = getMateCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBamForChr(INPUT, IntermediateBamType.MATE, dictionary.getSequence(i).getSequenceName()));
	    			SAMFileHeader mateHeader = header.clone();
	    			mateHeader.setSortOrder(SortOrder.unsorted);
	    			matecollection.add(SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(mateHeader), new SAMRecordMateCoordinateComparator(),
	    					// TODO: allocate buffers according to sequence lengths instead of equal space to every chr
	    					MAX_RECORDS_IN_RAM / dictionary.size(),
	    					TMP_DIR));
	    		}
				oeawriters[dictionary.size()] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBamForChr(INPUT, IntermediateBamType.OEA, "unmapped"));
	    	} else {
	    		// all writers map to the same one
	    		scwriters[0] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.SC));
				oeawriters[0] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.OEA));
				dpwriters[0] = getCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.DP));
				matewriters[0] = getMateCoordinateSortedWriter(header, FileNamingConvention.GetIntermediateBam(INPUT, IntermediateBamType.MATE));
				SAMFileHeader mateHeader = header.clone();
				mateHeader.setSortOrder(SortOrder.unsorted);
				matecollection.add(SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(mateHeader), new SAMRecordMateCoordinateComparator(), MAX_RECORDS_IN_RAM, TMP_DIR));
				for (int i = 1; i < dictionary.size() + 1; i++) {
					scwriters[i] = scwriters[0];
					oeawriters[i] = oeawriters[0];
					dpwriters[i] = dpwriters[0];
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
				if (SAMRecordSummary.isAlignmentSoftClipped(record)) {
					// soft clipped
					scwriters[offset].addAlignment(record);
				}
				if (record.getReadPairedFlag()) {
					if (!record.getProperPairFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
						// discordant pair
						dpwriters[offset].addAlignment(record);
						matecollection.get(record.getMateReferenceIndex()).add(record);
					}
					if (record.getReadUnmappedFlag() ^ record.getMateUnmappedFlag()) { // only 1 read of the pair is mapped
						// OEA
						oeawriters[offset].addAlignment(record);
						if (record.getReadUnmappedFlag()) {
							matecollection.get(record.getMateReferenceIndex()).add(record);
						}
					}
				}
				metrics.acceptRecord(record, null);
				progress.record(record);
			}
			reader.close();
			metrics.finish();
			final MetricsFile<InsertSizeMetrics, Integer> serialisedMetrics = getMetricsFile();
			metrics.addAllLevelsToFile(serialisedMetrics);
			serialisedMetrics.write(FileNamingConvention.GetMetrics(INPUT));
			for (int i = 0; i < dictionary.size(); i++) {
				scwriters[i].close();
				scwriters[i] = null;
				dpwriters[i].close();
				dpwriters[i] = null;
				oeawriters[i].close();
				oeawriters[i] = null;
			}
			if (PER_CHR) {
				for (int i = 0; i < dictionary.size(); i++) {
					matewriters[i].setProgressLogger(new ProgressLogger(log));
					writeToFile(matecollection.get(i), matewriters[i]);
					matewriters[i].close();
					matecollection.set(i, null);
				}
			} else {
				matewriters[0].setProgressLogger(new ProgressLogger(log));
				writeToFile(matecollection.get(0), matewriters[0]);
				matewriters[0].close();
				matecollection = null;
			}
    	} catch (IOException e) {
    		throw new RuntimeException(e);
    		
    	}
        return 0;
    }
	public static void main(String[] argv) {
        System.exit(new ExtractEvidence().instanceMain(argv));
    }
}
