package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.vcf.UntemplatedSequenceAnnotator;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

public class AnnotateUntemplatedSequence extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateUntemplatedSequence.class);
	@Argument(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF file to annotate")
    public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Annotated VCF file")
    public File OUTPUT;
	@Argument(doc="Overwrite existing annotation. Setting to 'true' will replace any existing BEALN annotations. " +
			"Setting to 'false' will generate annotations only for records without an existing BEALN annotation. ")
	public boolean OVERWRITE = false;
	@Argument(doc="Command line arguments to run external aligner. "
			+ "Aligner output must be written to stdout and the records MUST match the input fastq order."
			+ " The aligner must support using \"-\" as the input filename when reading from stdin."
			+ "Java argument formatting is used with %1$s being the fastq file to align, "
			+ "%2$s the reference genome, and %3$d the number of threads to use.", optional=true)
	public List<String> ALIGNER_COMMAND_LINE = new SoftClipsToSplitReads().ALIGNER_COMMAND_LINE;
	public static void main(String[] argv) {
        System.exit(new AnnotateUntemplatedSequence().instanceMain(argv));
    }
	@Override
	public int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		log.info("Annotating variant untemplated sequence in " + INPUT);
		try (UntemplatedSequenceAnnotator ann = new UntemplatedSequenceAnnotator(REFERENCE_SEQUENCE, INPUT, OVERWRITE, ALIGNER_COMMAND_LINE, WORKER_THREADS)) {
			saveVcf(INPUT, OUTPUT, ann);
			log.info("Annotated variants written to " + OUTPUT);
		} catch (IOException e) {
			log.error(e);
			throw new RuntimeException(e);
		}
		return 0;
	}
	protected void saveVcf(File input, File output, Iterator<VariantContext> calls) throws IOException {
		VCFHeader header;
		try (VCFFileReader vcfReader = new VCFFileReader(input, false)) {
			header = vcfReader.getFileHeader();
		}

		header.addMetaDataLine(VcfInfoAttributes.BREAKEND_ALIGNMENTS.infoHeader());
		File tmp = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(output) : output;
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
				.setOutputFile(tmp);
				//.setOption(Options.INDEX_ON_THE_FLY) // don't know if we're sorted or not so we can't index
		if (header.getSequenceDictionary() != null) {
			// annotation dictionary could be different to reference dictionary
			builder = builder.setReferenceDictionary(header.getSequenceDictionary());
		}
		final ProgressLogger writeProgress = new ProgressLogger(log);
		try (VariantContextWriter vcfWriter = builder.build()) {
			vcfWriter.writeHeader(header);
			while (calls.hasNext()) {
				VariantContext record = calls.next();
				vcfWriter.add(record);
				writeProgress.record(record.getContig(), record.getStart());
			}
		}
		if (tmp != output) {
			FileHelper.move(tmp, output, true);
		}
	}
}
