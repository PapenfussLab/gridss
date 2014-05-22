package au.edu.wehi.socrates;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

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
    @Option(doc="Reference used for alignment",
    		shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
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
	    	IOUtil.assertFileIsWritable(FileNamingConvention.getMetrics(INPUT));
	    	final SamReaderFactory samFactory = SamReaderFactory.makeDefault();
	    	final SamReader reader = samFactory.open(INPUT);
	    	final SAMFileHeader header = reader.getFileHeader();
	    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
	    	final ReferenceSequenceFileWalker referenceWalker = new ReferenceSequenceFileWalker(REFERENCE);
	    	final InsertSizeMetricsCollector metrics = RelevantMetrics.createCollector(header);
	    	final ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
	    	
	    	// Check we aligned to the reference supplied 
	    	if (reference.getSequenceDictionary() == null) {
	    		throw new IllegalArgumentException(String.format("Missing .dict file for reference %s. Please create using picard tools CreateSequenceDictionary", REFERENCE));
	    	}
	    	dictionary.assertSameDictionary(reference.getSequenceDictionary());
	    	
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
	    			writers[i] = factory.makeSAMOrBAMWriter(svHeader, true, FileNamingConvention.getSVBamForChr(INPUT, dictionary.getSequence(i).getSequenceName()));
	    			matewriters[i] = factory.makeSAMOrBAMWriter(mateHeader, false, FileNamingConvention.getMateBamForChr(INPUT, dictionary.getSequence(i).getSequenceName()));
	    			matecollection.add(SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(mateHeader), new SAMRecordMateCoordinateComparator(),
	    					// TODO: allocate buffers according to sequence lengths instead of equal space to every chr
	    					MAX_RECORDS_IN_RAM / dictionary.size(),
	    					TMP_DIR));
	    		}
				writers[dictionary.size()] = factory.makeSAMWriter(svHeader, true, FileNamingConvention.getSVBamForChr(INPUT, "unmapped"));
	    	} else {
	    		// all writers map to the same one
	    		writers[0] = factory.makeSAMOrBAMWriter(svHeader, true, FileNamingConvention.getSVBam(INPUT));
				matewriters[0] = factory.makeSAMOrBAMWriter(mateHeader, false, FileNamingConvention.getMateBam(INPUT));
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
				boolean sc = SAMRecordUtil.isAlignmentSoftClipped(record);
				boolean badpair = SAMRecordUtil.isPartOfNonReferenceReadPair(record);
				if (sc || badpair) {
					SAMRecordUtil.ensureNmTag(referenceWalker, record);
					writers[offset].addAlignment(record);
				}
				if (badpair) {
					if (!record.getMateUnmappedFlag()) {
						matecollection.get(record.getMateReferenceIndex()).add(record);
					}
				}
				metrics.acceptRecord(record, null);
				progress.record(record);
			}
			reader.close();
			log.info("Writing metrics");
			metrics.finish();
			RelevantMetrics.save(metrics, this.<InsertSizeMetrics, Integer>getMetricsFile(), FileNamingConvention.getMetrics(INPUT));
			for (int i = 0; i < dictionary.size(); i++) {
				if (writers[i] != null) {
					writers[i].close();
					writers[i] = null;
				}
			}
			log.info("Sorting sv mates");
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
