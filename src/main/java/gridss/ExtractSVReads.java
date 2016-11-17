package gridss;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import gridss.analysis.InsertSizeDistribution;
import gridss.filter.ClippedReadFilter;
import gridss.filter.IndelReadFilter;
import gridss.filter.OneEndAnchoredReadFilter;
import gridss.filter.ReadPairConcordanceFilter;
import gridss.filter.SplitReadFilter;
import gridss.filter.UnionAggregateFilter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Extracts reads and read pairs supporting putative structural variations.",
        usageShort = "Extracts reads and read pairs supporting putative structural variations."
)
public class ExtractSVReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(ExtractSVReads.class);
	private static final int ASYNC_BUFFERS = 2;
	private static final int ASYNC_BUFFER_SIZE = 300;
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file containing subset of input", optional=false)
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
    @Option(doc="Include read pairs in which only one of the read is aligned to the reference.", optional=true)
    public boolean SINGLE_MAPPED_PAIRED = true;
    @Option(doc="Include read pairs that align do not align in the expected orientation within the expected fragment size distribution.", optional=true)
    public boolean DISCORDANT_READ_PAIRS = true;
    @Option(doc="Method of calculating read pair concordance. Valid values are SAM_FLAG, PERCENTAGE, and FIXED", optional=true)
    public ReadPairConcordanceMethod READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.SAM_FLAG;
    @Option(doc="Minimum concordant read pair fragment size if using the FIXED method of calculation", optional=true)
    public int FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = 0;
    @Option(doc="Maximum concordant read pair fragment size if using the FIXED method of calculation", optional=true)
    public int FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = 0;
    @Option(doc = "Percent (0.0-1.0) of read pairs considered concorant if using the PERCENTAGE method of calculation.", optional=true)
    public Float READ_PAIR_CONCORDANT_PERCENT = 0.995f;
    @Option(doc="Picard tools insert size distribution metrics txt file. Required if using the PERCENTAGE read pair concordance calculation method.", optional=true)
    public File INSERT_SIZE_METRICS = null;
    @Option(doc="Include unmapped reads", optional=true)
    public boolean UNMAPPED_READS = true;
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	
    	// Read metrics file
    	InsertSizeDistribution isd = null;
    	if (INSERT_SIZE_METRICS != null) {
    		IOUtil.assertFileIsReadable(INSERT_SIZE_METRICS);
    		isd = InsertSizeDistribution.create(INSERT_SIZE_METRICS); 
    	}
    	ReadPairConcordanceCalculator rpcc = ReadPairConcordanceCalculator.create(
    			READ_PAIR_CONCORDANCE_METHOD,
    			FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE,
    			FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE,
    			READ_PAIR_CONCORDANT_PERCENT,
    			isd,
    			null);
    	
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
								SPLIT,
								SINGLE_MAPPED_PAIRED,
								DISCORDANT_READ_PAIRS, rpcc,
								UNMAPPED_READS);
    				}
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
    
	public static void extract(
			final Iterator<SAMRecord> rawit,
			final SAMFileWriter writer,
			final int minIndelSize,
			final int minClipLength,
			final boolean includeSplitReads,
			final boolean includeOEA,
			final boolean includeDP,
			final ReadPairConcordanceCalculator rpcc,
			final boolean includeUnmapped) throws IOException {
		ProgressLogger progress = new ProgressLogger(log);
		List<SamRecordFilter> filters = new ArrayList<>();
		filters.add(new IndelReadFilter(minIndelSize));
		filters.add(new ClippedReadFilter(minClipLength));
		if (includeSplitReads) filters.add(new SplitReadFilter());
		if (includeOEA) filters.add(new OneEndAnchoredReadFilter());
		if (includeDP) filters.add(new ReadPairConcordanceFilter(rpcc, false, true));
		if (includeUnmapped) filters.add(new AlignedFilter(false));
		
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
