package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.analysis.CollectStructuralVariantReadMetrics;
import gridss.cmdline.ProcessStructuralVariantReadsCommandLineProgram;
import gridss.filter.*;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

@CommandLineProgramProperties(
		summary = "Extracts reads and read pairs with a mapping to a region to extract.",
        oneLineSummary = "Extracts reads and read pairs with a mapping to a region to extract.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class ExtractReads extends CommandLineProgram {
	private static final Log log = Log.getInstance(ExtractReads.class);
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file")
	public File OUTPUT;
	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input file")
	public File INPUT;
	@Argument(shortName = "F", doc = "BED file containing regions to export")
	public File FILTER_BED;
	@Argument(doc = "Extract reads whose mate maps to an export region")
	public boolean EXTRACT_MATES = true;
	@Argument(doc = "Extract reads that have a chimeric alignment mapping to an export region")
	public boolean EXTRACT_SPLIT_READ = true;
	@Argument(doc = "Process the entire input file. A linear scan of the entire input file is" +
			"faster for Small files, files containing many export regions, " +
			"and sequencing data with high discordant read pair mapping rate.")
	public boolean LINEAR_SCAN = false;
	@Argument(doc = "Number of bases surrounding each export region to include in the index query. " +
			"This should be set to a value greater than the fragment size distribution.")
	public int REGION_PADDING_SIZE = 2000;

	public static void main(String[] argv) {
		System.exit(new ExtractReads().instanceMain(argv));
	}

	@Override
	protected int doWork() {
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsReadable(FILTER_BED);
		// create query intervals
			// Y/N lookup
			// Merge expanded intervals
		// Chunk into tasks
		// doPar
		{
			SamReader reader = SamReaderFactory.make().open(INPUT);
			// read padded interval
			// track out of region reads+position
			// unmapped should be in region (hopefully!)
		}
		// doPar
		// Extract out of region intervals/reads
		{
		}
	}
}
