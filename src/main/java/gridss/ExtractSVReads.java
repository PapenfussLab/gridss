package gridss;

import gridss.filter.ClippedReadFilter;
import gridss.filter.IndelReadFilter;
import gridss.filter.SplitReadFilter;
import gridss.filter.UnionAggregateFilter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;

/**
 * Extracts all reads and read pairs supporting any putative structural variation
 * @author Daniel Cameron
 *
 */
public class ExtractSVReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(ExtractSVReads.class);
	private static final int ASYNC_BUFFERS = 2;
	private static final int ASYNC_BUFFER_SIZE = 300;
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Input BEDPE", optional=false)
    public File OUTPUT;
    @Option(doc="Minimum indel size", optional=true)
    public int MIN_INDEL_SIZE = 1;
    @Option(doc="Minimum bases clipped", optional=true)
    public int MIN_CLIP_LENGTH = 1;
    @Option(doc="Include hard and soft clipped reads in output", optional=true)
    public boolean CLIPPED = true;
    @Option(doc="Include reads containing indels in output", optional=true)
    public boolean INDELS = true;
    @Option(doc="Include split reads in output", optional=true)
    public boolean SPLIT = true;
    // TODO: extract read pairs
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.make();
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	try {
    		try (SamReader reader = readerFactory.open(INPUT)) {
    			SAMFileHeader header = reader.getFileHeader();
    			try (SAMRecordIterator it = reader.iterator()) {
    				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, OUTPUT)) {
    					extract(it, 
    							writer,
    							INDELS ? MIN_INDEL_SIZE : Integer.MAX_VALUE,
    							CLIPPED ? MIN_CLIP_LENGTH : Integer.MAX_VALUE,
								SPLIT);
    				}
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	public static void extract(final Iterator<SAMRecord> rawit, final SAMFileWriter writer, final int minIndelSize, final int minClipLength, final boolean includeSplitReads) throws IOException {
		ProgressLogger progress = new ProgressLogger(log);
		List<SamRecordFilter> filters = new ArrayList<>();
		filters.add(new IndelReadFilter(minIndelSize));
		filters.add(new ClippedReadFilter(minClipLength));
		if (includeSplitReads) filters.add(new SplitReadFilter());
		UnionAggregateFilter filter = new UnionAggregateFilter(filters);
		try (CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(rawit, "raw", ASYNC_BUFFERS, ASYNC_BUFFER_SIZE)) {
			while (it.hasNext()) {
				SAMRecord r = it.next();
				progress.record(r);
				if (!filter.filterOut(r)) {
					writer.addAlignment(r);
				}
			}
		}
	}
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new ExtractSVReads().instanceMain(argv));
    }
}
