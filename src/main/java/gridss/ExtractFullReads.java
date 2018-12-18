package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.File;
import java.io.IOException;

@CommandLineProgramProperties(
		summary = "Extracts reads and read pairs with a mapping to a region to extract.",
        oneLineSummary = "Extracts reads and read pairs with a mapping to a region to extract.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class ExtractFullReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(ExtractFullReads.class);
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file")
	public File OUTPUT;
	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file")
	public File INPUT;
	@Argument(shortName = "B", doc = "BED file containing regions to export")
	public File REGION_BED;
	@Argument(doc = "Extract reads whose mate maps to an export region. " +
			"If the MC tag is not present, only the starting alignment position of the mate is considered. " +
			"When determining whether the mate maps to an export region " +
			"only the primary alignment of that mate is considered. Secondary " +
			"and supplementary alignments are ignored.")
	public boolean EXTRACT_MATES = true;
	@Argument(doc = "Extract all records for reads that have a chimeric alignment mapping to an export region")
	public boolean EXTRACT_SPLITS = true;
	@Argument(doc = "Process the entire input file. A linear scan of the entire input file is" +
			"faster for Small files, files containing many export regions, " +
			"and sequencing data with high discordant read pair mapping rate.")
	public boolean LINEAR_SCAN = false;
	@Argument(doc = "Number of bases surrounding each export region to include in the index query. " +
			"Setting this to a value slightly greater than the 99.99% fragment size distribution will reduce the number of random access" +
			"IO requests made. " +
			"This parameter is not used if a linear scan is performed.")
	public int REGION_PADDING_SIZE = 2000;

	@Argument(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
			shortName="THREADS")
	public int WORKER_THREADS = 1;

	public static void main(String[] argv) {
		System.exit(new ExtractFullReads().instanceMain(argv));
	}

	@Override
	protected int doWork() {
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsReadable(REGION_BED);
		try {
			SamReader reader = SamReaderFactory.make().open(INPUT);
			SAMFileHeader header = reader.getFileHeader();
			if (!LINEAR_SCAN) {
				if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
					log.error("INPUT is not coordinate sorted. "
							+ "Specify LINEAR_SCAN=true when processing files that are not coordinate sorted.");
					return -1;
				} else if (!reader.hasIndex()) {
					log.error("INPUT does not have an associated index. "
							+ "Specify LINEAR_SCAN=true or create an index file.");
					return -1;
				}
			}
			SAMSequenceDictionary dict = header.getSequenceDictionary();
			if (dict == null) {
				throw new IllegalArgumentException("Missing input file sequence dictionary. Ensure the correct @SQ headers exist.");
			}
			LinearGenomicCoordinate lgc = new PaddedLinearGenomicCoordinate(dict, GenomicProcessingContext.LINEAR_COORDINATE_CHROMOSOME_BUFFER, true);
			IntervalBed bed = new IntervalBed(dict, lgc, REGION_BED);

			FullReadExtractor extractor;
			if (LINEAR_SCAN) {
				extractor = new LinearScanFullReadExtractor(lgc, bed, EXTRACT_MATES, EXTRACT_SPLITS);
			} else {
				extractor = new IndexedLookupFullReadExtractor(lgc, bed, EXTRACT_MATES, EXTRACT_SPLITS, REGION_PADDING_SIZE);
			}
			extractor.extract(INPUT, OUTPUT, WORKER_THREADS);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return 0;
	}
}
