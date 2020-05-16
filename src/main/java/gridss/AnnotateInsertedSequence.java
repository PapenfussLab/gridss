package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.VariantContextRepeatMaskerAnnotator;
import au.edu.wehi.idsv.alignment.BwaStreamingAligner;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.vcf.InsertedSequenceAnnotator;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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
import java.util.Set;
import java.util.stream.Collectors;

public class AnnotateInsertedSequence extends ReferenceCommandLineProgram {
    private static final Log log = Log.getInstance(AnnotateInsertedSequence.class);
    @Argument(doc = "Number of worker threads to spawn. Defaults to number of cores available."
            + " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
            shortName = "THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "VCF file to annotate")
    public File INPUT;
    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Annotated VCF file")
    public File OUTPUT;
    @Argument(doc = "Command line arguments to run external aligner. "
            + "In-process bwa alignment is used if this value is null. "
            + "Aligner output must be written to stdout and the records MUST match the input fastq order."
            + " The aligner must support using \"-\" as the input filename when reading from stdin."
            + "Java argument formatting is used with %1$s being the fastq file to align, "
            + "%2$s the reference genome, and %3$d the number of threads to use.", optional = true)
    public List<String> ALIGNER_COMMAND_LINE = null;
    @Argument(doc = "Number of records to buffer when performing in-process or streaming alignment. Not applicable when performing external alignment.", optional = true)
    public int ALIGNER_BATCH_SIZE = MAX_RECORDS_IN_RAM;
    @Argument(doc = "Whether to align inserted sequences to REFERENCE_GENOME. Valid values are:" +
            "APPEND (Append alignments to REFERENCE_GENOME to the BEALN field), " +
            "REPLACE (Replace all BEALN fields)  (default)," +
            "ADD_MISSING (Add alignments to records missing a BEALN field, and" +
            "SKIP (do not align).", optional = true)
    public AlignmentStatus ALIGNMENT = AlignmentStatus.REPLACE;
    @Argument(doc = "Annotate inserted sequences with RepeatMasker annotations. Use bedops rmsk2bed to generate the bed file from the RepeatMasker .fa.out file.", optional = true)
    public File REPEAT_MASKER_BED = null;

    public static void main(String[] argv) {
        System.exit(new AnnotateInsertedSequence().instanceMain(argv));
    }

    public enum AlignmentStatus {
        /**
         * Append alignments to REFERENCE_GENOME to the BEALN field
         */
        APPEND,
        /**
         * Replace all BEALN fields
         */
        REPLACE,
        /**
         * Add alignments to REFERENCE_GENOME to the BEALN field
         */
        ADD_MISSING,
        /**
         * Do not perform additional alignments
         */
        SKIP;
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (REPEAT_MASKER_BED != null && !REPEAT_MASKER_BED.isFile()) {
            return new String[]{"REPEAT_MASKER_BED: file not found"};
        }
        return super.customCommandLineValidation();
    }

    @Override
    public int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        log.info("Annotating inserted sequences in " + INPUT);
        try {
            SAMSequenceDictionary dict = new IndexedFastaSequenceFile(REFERENCE_SEQUENCE).getSequenceDictionary();
            Iterator<VariantContext> it;
            if (ALIGNMENT != AlignmentStatus.SKIP) {
                StreamingAligner sa;
                if (ALIGNER_COMMAND_LINE == null || ALIGNER_COMMAND_LINE.size() == 0) {
                    log.info("Using in-process bwa alignment");
                    sa = new BwaStreamingAligner(REFERENCE_SEQUENCE, dict, WORKER_THREADS, ALIGNER_BATCH_SIZE * 150);
                } else {
                    log.info("Using external process alignment");
                    sa = new ExternalProcessStreamingAligner(SamReaderFactory.make(), ALIGNER_COMMAND_LINE, REFERENCE_SEQUENCE, WORKER_THREADS, dict);
                }
                InsertedSequenceAnnotator ann = new InsertedSequenceAnnotator(
                        INPUT,
                        sa,
                        ALIGNMENT == AlignmentStatus.REPLACE,
                        ALIGNMENT == AlignmentStatus.ADD_MISSING);
                it = ann;
            } else {
                VCFFileReader vcfReader = new VCFFileReader(INPUT, false);
                it = vcfReader.iterator();
            }
            if (REPEAT_MASKER_BED != null) {
                log.info("Loading RepeatMasker bed file from " + REPEAT_MASKER_BED);
                VariantContextRepeatMaskerAnnotator rma = new VariantContextRepeatMaskerAnnotator(REPEAT_MASKER_BED);
                Set<String> commonContigs = Sets.intersection(Sets.newHashSet(rma.getRepeatMaskerContigs()), dict.getSequences().stream().map(s -> s.getContig()).collect(Collectors.toSet()));
                if (commonContigs.size() < Math.min(rma.getRepeatMaskerContigs().size(), dict.size()) * 0.5) {
                    log.warn(String.format("Only %d chromosomes in common between REFERENCE_SEQUENCE and REPEAT_MASKER_BED. Are you sure your chromosome names match?", commonContigs.size()));
                }
                VariantContextRepeatMaskerAnnotator finalRma = rma;
                it = Iterators.transform(it, (vc) -> finalRma.apply(vc));
            }
            saveVcf(INPUT, OUTPUT, it);
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
        if (REPEAT_MASKER_BED != null) {
            header.addMetaDataLine(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_OVERLAP.infoHeader());
            header.addMetaDataLine(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_TYPE.infoHeader());
            header.addMetaDataLine(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_REPEAT_CLASS.infoHeader());
            header.addMetaDataLine(VcfInfoAttributes.INSERTED_SEQUENCE_REPEATMASKER_ORIENTATION.infoHeader());
        }
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
