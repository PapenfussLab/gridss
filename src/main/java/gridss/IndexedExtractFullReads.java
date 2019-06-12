package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

@CommandLineProgramProperties(
		summary = "Extracts reads and read pairs with a mapping to a region to extract from an indexed bam file. " +
				" Uses the bam file index. Note that since the read mate and split reads are also extracted, " +
				" many more locations that just those supplied in the bed file will be queried. " +
				" When querying small regions, an indexed scan may be faster. For extracting large " +
				" regions of the genome, the linear scanning of ExtractFullReads will likely be faster.",
        oneLineSummary = "Extracts reads and read pairs with a mapping to a region to extract from an indexed bam file.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class IndexedExtractFullReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(IndexedExtractFullReads.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Index bam file.")
	public File INPUT;
	@Argument(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file location.")
	public File OUTPUT;
	@Argument(shortName = "B", doc = "BED file containing regions to export")
	public File REGION_BED = new ExtractFullReads().REGION_BED;
	@Argument(doc = "Extract reads whose mate maps to an export region. " +
			"If the MC tag is not present, only the starting alignment position of the mate is considered. " +
			"When determining whether the mate maps to an export region " +
			"only the primary alignment of that mate is considered. Secondary " +
			"and supplementary alignments are ignored.")
	public boolean EXTRACT_MATES = new ExtractFullReads().EXTRACT_MATES;
	@Argument(doc = "Extract all records for reads that have a chimeric alignment mapping to an export region")
	public boolean EXTRACT_SPLITS = new ExtractFullReads().EXTRACT_SPLITS;
	@Argument(doc = "Number of additional bases surrounding each export region to include in the index query. ")
	public int REGION_PADDING_SIZE = new ExtractFullReads().REGION_PADDING_SIZE;
	@Argument(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
			shortName="THREADS")
	public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();

	public static void main(String[] argv) {
		System.exit(new IndexedExtractFullReads().instanceMain(argv));
	}

	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
		java.util.Locale.setDefault(Locale.ROOT);
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsReadable(REGION_BED);
		try {
			SAMSequenceDictionary dict = null;
			try (SamReader reader = SamReaderFactory.makeDefault().open(INPUT)) {
				if (!reader.hasIndex()) {
					log.error("Missing BAM index for " + INPUT.getName());
					return -1;
				}
				dict = reader.getFileHeader().getSequenceDictionary();
			}
			if (dict == null) {
				log.error("Missing input file sequence dictionary. Ensure the correct @SQ headers exist.");
				return -1;
			}
			if (dict == null) {
				throw new IllegalArgumentException("Missing input file sequence dictionary. Ensure the correct @SQ headers exist.");
			}
			LinearGenomicCoordinate lgc = new PaddedLinearGenomicCoordinate(dict, GenomicProcessingContext.LINEAR_COORDINATE_CHROMOSOME_BUFFER, true);
			IntervalBed bed = new IntervalBed(lgc, REGION_BED);
			bed = bed.expandIntervals(REGION_PADDING_SIZE, REGION_PADDING_SIZE);
			ReadExtractor extractor = new IndexedReadExtractor(lgc, bed, EXTRACT_MATES, EXTRACT_SPLITS);
			extractor.extract(INPUT, OUTPUT, WORKER_THREADS);
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
}
