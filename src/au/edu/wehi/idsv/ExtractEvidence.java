package au.edu.wehi.idsv;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Extracts reads supporting structural variation into intermediate files.
 * By sorting DP and OEA read pairs by the coordinate of the mapped mate,
 * putative directed breakpoints can be assembled by downstream processes
 * in a single pass of the intermediate files.
 * @author Daniel Cameron
 *
 */
public class ExtractEvidence implements Closeable {
	public static int MAX_RECORDS_IN_RAM = 500000; // match htsjdk default
	private static final Log log = Log.getInstance(ExtractEvidence.class);
	private final ProcessingContext processContext;
	private final File input;
	private SamReader reader = null;
	private ReferenceSequenceFileWalker referenceWalker = null;
	private List<SAMFileWriter> writers = Lists.newArrayList();
	private List<SAMFileWriter> matewriters = Lists.newArrayList();
	private List<SortingCollection<SAMRecord>> mateRecordBuffer = Lists.newArrayList();
	/**
	 * Files to delete if there is an error during processing.
	 * Deleting output files if there is an error prevents downstream issues
	 * with partially written intermediate files 
	 */
	private Set<File> deleteOnError = Sets.newHashSet();
	public void close() {
		tryClose(reader);
		reader = null;
		tryClose(referenceWalker);
		referenceWalker = null;
		for (SAMFileWriter w : writers) {
			tryClose(w);
		}
    	writers = null;
    	for (SAMFileWriter w : matewriters) {
			tryClose(w);
		}
    	matewriters = null;
	}
	private void tryClose(Closeable toClose) {
		try {
			if (toClose != null) toClose.close();
    	} catch (IOException e) {
    		log.equals(e);
    	}
	}
	/**
	 * Deletes all output files
	 */
	private void deleteOutput() {
		close(); // close any file handles that are still around
		for (File f : deleteOnError) {
			f = new File(f.getAbsolutePath());
			if (f.exists()) {
				try {
					if (!f.delete()) {
						log.error("Unable to delete %s during output rollback", f);
					}
				} catch (Exception e) {
					log.error(e);
				}
			}
		}
	}
    private void writeToFile(SortingCollection<SAMRecord> sc, SAMFileWriter writer) {
    	sc.doneAdding();
    	final CloseableIterator<SAMRecord> it = sc.iterator();
		while (it.hasNext()) {
			writer.addAlignment(it.next());
		}
    }
    public ExtractEvidence(ProcessingContext processContext, File input) {
    	this.processContext = processContext;
    	this.input = input;
    }
	public void extractEvidence() {
    	try {
	    	IOUtil.assertFileIsWritable(processContext.getFileSystemContext().getMetrics(input));
	    	
	    	reader = processContext.getSamReaderFactory().open(input);
	    	final SAMFileHeader header = reader.getFileHeader();
	    	final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
	    	referenceWalker = new ReferenceSequenceFileWalker(processContext.getReferenceFile());
	    	final InsertSizeMetricsCollector metrics = RelevantMetrics.createCollector(header);
	    	
	    	dictionary.assertSameDictionary(processContext.getReference().getSequenceDictionary());
	    	
	    	createSamOutputWriters(header, dictionary);
	    	
	    	processInputRecords(dictionary, metrics);
	    	
	    	writeMetrics(metrics);
	    	
			flush();
    	} catch (Exception e) {
    		deleteOutput();
    		throw new RuntimeException(e);
    	} finally {
    		close();
    	}
    }
	private void writeMetrics(InsertSizeMetricsCollector metrics) throws IOException {
		log.info("Writing metrics");
		metrics.finish();
		RelevantMetrics.save(metrics, processContext.<InsertSizeMetrics, Integer>createMetricsFile(), processContext.getFileSystemContext().getMetrics(input));
	}
	private void flush() throws IOException {
		reader.close();
		reader = null;
		for (SAMFileWriter w : writers) {
			w.close();
		}
		writers.clear();
		log.info("Sorting sv mates");
		for (int i = 0; i < mateRecordBuffer.size(); i++) {
			// TODO: multi-thread this?
			SortingCollection<SAMRecord> buffer = mateRecordBuffer.get(i);
			SAMFileWriter w = matewriters.get(i);
			w.setProgressLogger(new ProgressLogger(log));
			writeToFile(buffer, w);
			w.close();
			matewriters.set(i, null);
			mateRecordBuffer.set(i,  null);
		}
		mateRecordBuffer.clear();
		matewriters.clear();
	}
	private void processInputRecords(final SAMSequenceDictionary dictionary, final InsertSizeMetricsCollector metrics) {
		// output all to a single file, or one per chr + one for unmapped 
		assert(matewriters.size() == mateRecordBuffer.size());
		assert(writers.size() == dictionary.getSequences().size() + 1 || writers.size() == 1);
		assert(matewriters.size() == dictionary.getSequences().size() || matewriters.size() == 1);
		
		// Traverse the input file
		final ProgressLogger progress = new ProgressLogger(log);
		final SAMRecordIterator iter = reader.iterator().assertSorted(SortOrder.coordinate);
		while (iter.hasNext()) {
			SAMRecord record = iter.next();
			int offset = record.getReadUnmappedFlag() ? dictionary.size() : record.getReferenceIndex();
			boolean sc = SAMRecordUtil.isSoftClipLengthAtLeast(record, processContext.getSoftClipParameters().minLength);
			boolean badpair = SAMRecordUtil.isPartOfNonReferenceReadPair(record);
			if (sc || badpair) {
				SAMRecordUtil.ensureNmTag(referenceWalker, record);
				writers.get(offset % writers.size()).addAlignment(record);
			}
			if (badpair) {
				if (!record.getMateUnmappedFlag()) {
					mateRecordBuffer.get(record.getMateReferenceIndex() % mateRecordBuffer.size()).add(record);
				}
			}
			metrics.acceptRecord(record, null);
			progress.record(record);
		}
	}
	private void createSamOutputWriters(final SAMFileHeader header, final SAMSequenceDictionary dictionary) throws IOException {
		final SAMFileHeader svHeader = header.clone();
		svHeader.setSortOrder(SortOrder.coordinate);
		final SAMFileHeader mateHeader = header.clone();
		mateHeader.setSortOrder(SortOrder.unsorted);
		if (processContext.shouldProcessPerChromosome()) {
			createSamOutputWritersPerChromosome(dictionary, svHeader, mateHeader);
		} else {
			createSamOutputWriterPerGenome(svHeader, mateHeader);
		}
	}
	private void createSamOutputWriterPerGenome(final SAMFileHeader svHeader, final SAMFileHeader mateHeader) throws IOException {
		// all writers map to the same one
		writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, track(processContext.getFileSystemContext().getSVBam(input))));
		matewriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(mateHeader, false, track(processContext.getFileSystemContext().getMateBam(input))));
		mateRecordBuffer.add(SortingCollection.newInstance(
				SAMRecord.class,
				new BAMRecordCodec(mateHeader),
				new SAMRecordMateCoordinateComparator(),
				processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
				processContext.getFileSystemContext().getTemporaryDirectory()));
	}
	private void createSamOutputWritersPerChromosome(
			final SAMSequenceDictionary dictionary,
			final SAMFileHeader svHeader, final SAMFileHeader mateHeader)
			throws IOException {
		for (SAMSequenceRecord seq : dictionary.getSequences()) {
			writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, track(processContext.getFileSystemContext().getSVBamForChr(input, seq.getSequenceName()))));
			matewriters.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(mateHeader, false, track(processContext.getFileSystemContext().getMateBamForChr(input, seq.getSequenceName()))));
			mateRecordBuffer.add(SortingCollection.newInstance(
					SAMRecord.class,
					new BAMRecordCodec(mateHeader),
					new SAMRecordMateCoordinateComparator(),
					processContext.getFileSystemContext().getMaxBufferedRecordsPerFile(),
					processContext.getFileSystemContext().getTemporaryDirectory()));
		}
		writers.add(processContext.getSamReaderWriterFactory().makeSAMOrBAMWriter(svHeader, true, track(processContext.getFileSystemContext().getSVBamForChr(input, "unmapped"))));
	}
	/**
	 * Tracks the given output file, deleting it if there is an error. 
	 * @param output file to track
	 * @return tracked file
	 */
	public File track(File output) {
		deleteOnError.add(output);
		return output;
	}
}
