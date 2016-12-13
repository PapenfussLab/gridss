package gridss;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import com.google.common.base.Strings;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProgressLoggingSAMRecordIterator;
import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.analysis.InsertSizeDistribution;
import gridss.filter.ClippedReadFilter;
import gridss.filter.FixedFilter;
import gridss.filter.IndelReadFilter;
import gridss.filter.OneEndAnchoredReadFilter;
import gridss.filter.ReadPairConcordanceFilter;
import gridss.filter.SplitReadFilter;
import gridss.filter.UnionAggregateFilter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
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
        usage = "Extracts reads and read pairs supporting putative structural variations. "
        		+ "If the input file is queryname sorted, a multi-mapping aware extraction is performed "
        		+ "and reads/read pairs are only extracted when all alignments are consistent with the "
        		+ "presence of of a structural variant.",
        usageShort = "Extracts reads and read pairs supporting putative structural variations."
)
public class ExtractSVReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(ExtractSVReads.class);
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file. "
    		+ "If multiple mapping locations are reported for each read, these reads must be grouped together "
    		+ "(eg name sorted).", optional=false)
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
    		if (!INSERT_SIZE_METRICS.exists()) {
    			log.warn("Missing " + INSERT_SIZE_METRICS);
    		} else {
    			isd = InsertSizeDistribution.create(INSERT_SIZE_METRICS);
    		}
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
    			final SAMFileHeader header = reader.getFileHeader();
    			if (header.getSortOrder() != SortOrder.queryname) {
    				log.info("Not considering multiple read alignments as the input file is not queryname sorted.");
    			}
    			//final SAMSequenceDictionary dict = header.getSequenceDictionary();
    			//final LinearGenomicCoordinate linear = new PaddedLinearGenomicCoordinate(dict, 0);
    			//final SequentialCoverageThreshold coverageThresholdBlacklist = new SequentialCoverageThreshold(dict, linear, MAX_COVERAGE + 1);    			
    			try (SAMRecordIterator it = reader.iterator()) {
    				final File tmpoutput = FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ExtractSVReads.");
    				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
    					long count = extract(it, 
    							writer,
    							INDELS ? MIN_INDEL_SIZE : Integer.MAX_VALUE,
    							CLIPPED ? MIN_CLIP_LENGTH : Integer.MAX_VALUE,
								SPLIT,
								SINGLE_MAPPED_PAIRED,
								DISCORDANT_READ_PAIRS, rpcc,
								UNMAPPED_READS,
								INPUT.getName() + "-");
    					log.info(String.format("Extracted %d reads from %s", count, INPUT));
    				}
    				FileHelper.move(tmpoutput, OUTPUT, true);
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
    
	public static long extract(
			final Iterator<SAMRecord> rawit,
			final SAMFileWriter writer,
			final int minIndelSize,
			final int minClipLength,
			final boolean includeSplitReads,
			final boolean includeOEA,
			final boolean includeDP,
			final ReadPairConcordanceCalculator rpcc,
			final boolean includeUnmapped,
			final String threadPrefix) throws IOException {
		long count = 0;
		ProgressLogger progress = new ProgressLogger(log);
		List<SamRecordFilter> readfilters = new ArrayList<>();
		readfilters.add(new IndelReadFilter(minIndelSize));
		readfilters.add(new ClippedReadFilter(minClipLength));
		if (includeSplitReads) readfilters.add(new SplitReadFilter());
		if (includeUnmapped) readfilters.add(new AlignedFilter(false));
		SamRecordFilter readfilter = new UnionAggregateFilter(readfilters);
		
		List<SamRecordFilter> pairfilters = new ArrayList<>();
		if (includeOEA) pairfilters.add(new OneEndAnchoredReadFilter());
		if (includeDP) pairfilters.add(new ReadPairConcordanceFilter(rpcc, false, true));
		SamRecordFilter pairfilter = new UnionAggregateFilter(pairfilters);
		if (!includeOEA && !includeDP) {
			pairfilter = new FixedFilter(true);
		}

		try (CloseableIterator<SAMRecord> asyncit = new AsyncBufferedIterator<SAMRecord>(rawit, threadPrefix + "extract")) {
			PeekingIterator<SAMRecord> pit = Iterators.peekingIterator(new ProgressLoggingSAMRecordIterator(asyncit, progress));
			ArrayList<SAMRecord> records = new ArrayList<>();
			while (pit.hasNext()) {
				String readname = pit.peek().getReadName();
				records.clear();
				if (Strings.isNullOrEmpty(readname)) {
					// no read name so we just have to treat it as a single record
					records.add(pit.next());
				} else {
					while (readname != null && pit.hasNext() && readname.equals(pit.peek().getReadName())) {
						records.add(pit.next());
					}
				}
				boolean hasConsistentReadPair = hasReadPairingConsistentWithReference(rpcc, records);
				boolean[] hasConsistentReadAlignment = hasReadAlignmentConsistentWithReference(records);
				for (SAMRecord r : records) {
					if ((!readfilter.filterOut(r) && !hasConsistentReadAlignment[SAMRecordUtil.getSegmentIndex(r)])
							|| (!pairfilter.filterOut(r) && !hasConsistentReadPair)) {
						writer.addAlignment(r);
						count++;
					} else {
						// ignore remaining reads
					}
				}
			}
		}
		return count;
	}
	public static boolean[] hasReadAlignmentConsistentWithReference(List<SAMRecord> records) {
		boolean[] consistent = new boolean[2];
		for (SAMRecord r : records) {
			int segmentIndex = SAMRecordUtil.getSegmentIndex(r);
			if (consistent.length <= segmentIndex) {
				consistent = Arrays.copyOf(consistent, segmentIndex + 1);
			}
			if (consistent[segmentIndex]) continue; // no need to check this segment any further
			if (isFullyMapped(r)) {
				consistent[segmentIndex] = true;
			}
		}
		return consistent;
	}
	private static boolean isFullyMapped(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return false;
		Cigar cigar = r.getCigar();
		if (cigar == null) return false;
		for (CigarElement ce : cigar.getCigarElements()) {
			switch (ce.getOperator()) {
				case M:
				case EQ:
				case X:
				case P:
					break;
				default:
					return false;
			}
		}
		return true;
	}
	public static boolean hasReadPairingConsistentWithReference(ReadPairConcordanceCalculator rpcc, List<SAMRecord> records) {
		if (rpcc == null) return false;
		for (SAMRecord r1 : records) {
			if (r1.getReadUnmappedFlag() || !r1.getReadPairedFlag() || !r1.getFirstOfPairFlag()) continue;
			for (SAMRecord r2 : records) {
				if (r2.getReadUnmappedFlag() || !r2.getReadPairedFlag() || !r2.getSecondOfPairFlag()) continue;
				if (rpcc.isConcordant(r1, r2)) {
					return true;
				}
			}
		}
		return false;
	}

	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	@Override
	protected String[] customCommandLineValidation() {
		if (READ_PAIR_CONCORDANCE_METHOD == ReadPairConcordanceMethod.PERCENTAGE && INSERT_SIZE_METRICS == null) {
			return new String[] { "INSERT_SIZE_METRICS is required when using percentage based read pair concordance" };
		}
		return super.customCommandLineValidation();
	}
	public static void main(String[] argv) {
        System.exit(new ExtractSVReads().instanceMain(argv));
    }
}
