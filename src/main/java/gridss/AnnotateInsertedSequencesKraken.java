package gridss;

import au.edu.wehi.idsv.kraken.KrakenAnnotateVcf;
import au.edu.wehi.idsv.kraken.KrakenParser;
import au.edu.wehi.idsv.vcf.GridssVcfConstants;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;

import java.io.*;

@CommandLineProgramProperties(
        summary = "Annotate single breakends and breakend inserted sequences with a NCBI taxonomy ID. " +
                "Requires a pipe both to and from kraken2. Use mkfifo(1) to construct.",
        oneLineSummary = "Annotate VCFs and extract reads based on kraken2 classification",
        programGroup=gridss.cmdline.programgroups.DataConversion.class
)
public class AnnotateInsertedSequencesKraken extends CommandLineProgram {
    private static final Log log = Log.getInstance(AnnotateInsertedSequencesKraken.class);
    @Argument(doc="Output to be fed as the input to Kraken2. Should be a pipe.")
    public File FASTQ_PIPE_TO_KRAKEN;
    @Argument(doc="Kraken2 output to be read by this program.")
    public File OUTPUT_PIPE_FROM_KRAKEN;
    @Argument(doc="Minimum length of sequence to perform classification on.", optional=true)
    public int MIN_SEQUENCE_LENGTH = 20;
    @Argument(doc="GRIDSS VCF of calls generated from the input file.")
    public File INPUT_VCF;
    @Argument(doc="GRIDSS VCF annotated with kraken taxonomic identifiers")
    public File OUTPUT_VCF;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT_VCF);
        IOUtil.assertFileIsWritable(OUTPUT_VCF);
        try (KrakenParser parser = new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(OUTPUT_PIPE_FROM_KRAKEN))))) {
            try (BasicFastqWriter fqwriter = new BasicFastqWriter(new PrintStream(new FileOutputStream(FASTQ_PIPE_TO_KRAKEN, true)))) {
                log.info("Annotating NCBI Taxonomy IDs for " + INPUT_VCF);
                annotateVcf(parser, fqwriter, INPUT_VCF, OUTPUT_VCF, MIN_SEQUENCE_LENGTH);
            }
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    private static void annotateVcf(KrakenParser parser, BasicFastqWriter fqwriter, File inputVcf, File outputVcf, int minSequenceLength) {
        try (VCFFileReader vcfReader = new VCFFileReader(inputVcf, false)) {
            VCFHeader header = vcfReader.getFileHeader();
            SAMSequenceDictionary dict = header.getSequenceDictionary();
            CloseableIterator<VariantContext> it = vcfReader.iterator();
            VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                    .setOutputFile(outputVcf)
                    .setReferenceDictionary(dict);
            try (VariantContextWriter vcfWriter = builder.build()) {
                GridssVcfConstants.addHeaders(header);
                vcfWriter.writeHeader(header);
                KrakenAnnotateVcf kav = new KrakenAnnotateVcf(parser, fqwriter, it, dict, minSequenceLength);
                while (kav.hasNext()) {
                    VariantContext vc = kav.next();
                    vcfWriter.add(vc);
                }
            }
        }
    }

    public static void main(String[] argv) {
        System.exit(new AnnotateInsertedSequencesKraken().instanceMain(argv));
    }
}
