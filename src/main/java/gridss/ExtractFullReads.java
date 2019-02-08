package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.util.FileHelper;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.analysis.SinglePassSamProgram;

import java.io.File;
import java.io.IOException;

@CommandLineProgramProperties(
		summary = "Extracts reads and read pairs with a mapping to a region to extract.",
        oneLineSummary = "Extracts reads and read pairs with a mapping to a region to extract.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class ExtractFullReads extends SinglePassSamProgram {
	private static final Log log = Log.getInstance(ExtractFullReads.class);
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
	@Argument(doc = "Number of bases surrounding each export region to include in the index query. ")
	public int REGION_PADDING_SIZE = 2000;

	private File tmpOut;
	private SAMFileWriter writer;
	private FullReadExtractor extractor;

	public static void main(String[] argv) {
		System.exit(new ExtractFullReads().instanceMain(argv));
	}

	@Override
	protected void setup(SAMFileHeader header, File samFile) {
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsReadable(REGION_BED);
		SAMSequenceDictionary dict = header.getSequenceDictionary();
		if (dict == null) {
			throw new IllegalArgumentException("Missing input file sequence dictionary. Ensure the correct @SQ headers exist.");
		}
		LinearGenomicCoordinate lgc = new PaddedLinearGenomicCoordinate(dict, GenomicProcessingContext.LINEAR_COORDINATE_CHROMOSOME_BUFFER, true);
		tmpOut = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(OUTPUT) : OUTPUT;
		IntervalBed bed;
		try {
			bed = new IntervalBed(lgc, REGION_BED);
			bed = bed.expandIntervals(REGION_PADDING_SIZE, REGION_PADDING_SIZE);
		} catch (IOException e) {
			throw new RuntimeIOException(e);
		}
		extractor = new FullReadExtractor(lgc, bed, EXTRACT_MATES, EXTRACT_SPLITS);
		writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, tmpOut);
	}

	@Override
	protected void acceptRead(SAMRecord rec, ReferenceSequence ref) {
		if (extractor.shouldExtract(rec)) {
			writer.addAlignment(rec);
		}
	}

	@Override
	protected void finish() {
		writer.close();
		if (tmpOut != OUTPUT) {
			try {
				FileHelper.move(tmpOut, OUTPUT, true);
			} catch (IOException e) {
				throw new RuntimeIOException(e);
			}
		}
	}
}
