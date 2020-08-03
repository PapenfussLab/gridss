package gridss;

import au.edu.wehi.idsv.SAMRecordCoordinateOnlyComparator;
import au.edu.wehi.idsv.kraken.KrakenAnnotateVcf;
import au.edu.wehi.idsv.kraken.KrakenExtractReadInTaxonomy;
import au.edu.wehi.idsv.kraken.KrakenParser;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import au.edu.wehi.idsv.vcf.GridssVcfConstants;
import com.google.common.collect.Lists;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import joptsimple.internal.Strings;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Program to extract reads associated with specific NCBI taxonomy IDs that do not map to the reference genome." +
                "and/or annotate single breakends and breakend inserted sequences with a NCBI taxonomy ID. " +
                "Requires a pipe both to and from kraken2. Use mkfifo(1) to construct.",
        oneLineSummary = "Annotate VCFs and extract reads based on kraken2 classification",
        programGroup=gridss.cmdline.programgroups.DataConversion.class
)
public class ExtractFromTaxonomy extends CommandLineProgram {
    private static final int NCBI_VIRUS_TAXID = 10239;
    private static final Log log = Log.getInstance(ExtractFromTaxonomy.class);
    @Argument(doc="Output to be fed as the input to Kraken2. Should be a pipe.")
    public File FASTQ_PIPE_TO_KRAKEN;
    @Argument(doc="Kraken2 output to be read by this program.")
    public File OUTPUT_PIPE_FROM_KRAKEN;
    @Argument(doc="Input file to extract reads in which unmapped sequence maps to a NCBI taxonomy ID of interest.")
    public File EXTRACTION_INPUT_BAM;
    @Argument(doc="NCBI Taxonomy IDs to extract. All taxonomic entries under these IDs are also extracted. Defaults to all viruses.", optional = true)
    public List<Integer> EXTRACTION_TAXONOMY_IDS = Lists.newArrayList(NCBI_VIRUS_TAXID);
    @Argument(doc="NCBI taxonomy nodes.dmp. Download and extract from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
    public File NCBI_NODES_DMP;
    @Argument(doc="File to write read names of extracted reads.", optional=true)
    public File EXTRACTION_READ_LIST_FILE;
    @Argument(doc="File to extract unpaired reads to. Requires second pass over input file.", optional=true)
    public File EXTRACTION_OUTPUT_FQ;
    @Argument(doc="File to extract first read in pair to. Requires second pass over input file.", optional=true)
    public File EXTRACTION_OUTPUT_FQ1;
    @Argument(doc="File to extract second read in pair to. Requires second pass over input file.", optional=true)
    public File EXTRACTION_OUTPUT_FQ2;
    @Argument(doc="Minimum length of sequence to perform classification on.", optional=true)
    public int MIN_SEQUENCE_LENGTH = 20;
    @Argument(doc="Extract fragments mapping to the reference genome", optional=true)
    public boolean INCLUDE_MAPPED_FRAGMENTS_READS = false;
    @Argument(doc="Check taxonomy of soft clipped sequence when considering candidates for extraction.", optional=true)
    public boolean CONSIDER_SOFT_CLIPPED_BASES = true;

    @Override
    protected String[] customCommandLineValidation() {
        if (NCBI_NODES_DMP == null) {
            return new String[] { "NCBI_NODES_DMP required if EXTRACTION_INPUT_BAM is specified." };
        }
        if (EXTRACTION_TAXONOMY_IDS == null || EXTRACTION_TAXONOMY_IDS.isEmpty()) {
            return new String[] { "EXTRACTION_TAXONOMY_IDS required." };
        }
        if (EXTRACTION_READ_LIST_FILE == null && EXTRACTION_OUTPUT_FQ == null) {
            return new String[] { "EXTRACTION_READ_LIST_FILE or EXTRACTION_OUTPUT_FQ must be specified." };
        }
        if (!((EXTRACTION_OUTPUT_FQ == null && EXTRACTION_OUTPUT_FQ1 == null && EXTRACTION_OUTPUT_FQ2 == null) ||
                (EXTRACTION_OUTPUT_FQ != null && EXTRACTION_OUTPUT_FQ1 != null && EXTRACTION_OUTPUT_FQ2 != null))) {
            return new String[] { "EXTRACTION_OUTPUT_FQ, EXTRACTION_OUTPUT_FQ1, and EXTRACTION_OUTPUT_FQ2 must be specified together." };
        }
        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(EXTRACTION_INPUT_BAM);
        for (File outFile : Lists.newArrayList(
                EXTRACTION_INPUT_BAM,
                EXTRACTION_READ_LIST_FILE,
                EXTRACTION_OUTPUT_FQ,
                EXTRACTION_OUTPUT_FQ1,
                EXTRACTION_OUTPUT_FQ2)) {
            if (outFile != null) {
                IOUtil.assertFileIsWritable(outFile);
            }
        }
        try (KrakenParser parser = new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(OUTPUT_PIPE_FROM_KRAKEN))))) {
            try (BasicFastqWriter fqWriter = new BasicFastqWriter(new PrintStream(new FileOutputStream(FASTQ_PIPE_TO_KRAKEN, true)))) {
                log.info("Loading NCBI taxonomy from ", NCBI_NODES_DMP);
                boolean[] taxIdLookup = TaxonomyHelper.createInclusionLookup(EXTRACTION_TAXONOMY_IDS, TaxonomyHelper.parseMinimal(NCBI_NODES_DMP));
                log.info("Performing taxonomy lookup on ", EXTRACTION_INPUT_BAM);
                Set<String> readNames = new HashSet<>();
                try (SamReader reader = SamReaderFactory.makeDefault().open(EXTRACTION_INPUT_BAM)) {
                    KrakenExtractReadInTaxonomy kerit = new KrakenExtractReadInTaxonomy(parser, fqWriter, reader.iterator(), MIN_SEQUENCE_LENGTH, CONSIDER_SOFT_CLIPPED_BASES, INCLUDE_MAPPED_FRAGMENTS_READS, taxIdLookup);
                    while (kerit.hasNext()) {
                        String readName = kerit.next();
                        if (!Strings.isNullOrEmpty(readName)) {
                            readNames.add(readName);
                        }
                    }
                }
                if (EXTRACTION_READ_LIST_FILE != null) {
                    log.info("Writing matching read names to ", EXTRACTION_READ_LIST_FILE);
                    Files.write(EXTRACTION_READ_LIST_FILE.toPath(), readNames, StandardCharsets.UTF_8);
                }
                if (EXTRACTION_OUTPUT_FQ != null) {
                    log.info("Exporting matching reads to fastq.");
                    toFastq(EXTRACTION_INPUT_BAM, EXTRACTION_OUTPUT_FQ1, EXTRACTION_OUTPUT_FQ2, EXTRACTION_OUTPUT_FQ, readNames);
                }
            }
        } catch (FileNotFoundException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }
    private static FastqRecord samToFastq(SAMRecord r) {
        byte[] bases = r.getReadBases().clone();
        byte[] quals = r.getBaseQualities().clone();
        if (r.getReadNegativeStrandFlag()) {
            SequenceUtil.reverseComplement(bases);
            SequenceUtil.reverseQualities(quals);
        }
        return new FastqRecord(
                r.getReadName(),
                new String(bases, StandardCharsets.UTF_8),
                null,
                SAMUtils.phredToFastq(quals));
    }

    private void toFastq(
            File extraction_input_bam,
            File extraction_output_fq1,
            File extraction_output_fq2,
            File extraction_output_unpaired,
            Set<String> readNames) throws IOException {
        Map<String, SAMRecord> lookup = new HashMap<>();
        FastqWriterFactory factory = new FastqWriterFactory();
        try (FastqWriter fq1 = factory.newWriter(extraction_output_fq1)) {
            try (FastqWriter fq2 = factory.newWriter(extraction_output_fq2)) {
                try (FastqWriter fq = factory.newWriter(extraction_output_unpaired)) {
                    try (SamReader reader = SamReaderFactory.makeDefault().open(extraction_input_bam)) {
                        try (SAMRecordIterator it = reader.iterator()) {
                            while (it.hasNext()) {
                                if (readNames.isEmpty()) {
                                    log.debug("Found all reads. Stopping input file traversal.");
                                    break;
                                }
                                SAMRecord r = it.next();
                                String name = r.getReadName();
                                if (!r.getSupplementaryAlignmentFlag() && !r.isSecondaryAlignment() && readNames.contains(name)) {
                                    if (!r.getReadPairedFlag()) {
                                        fq.write(samToFastq(r));
                                        readNames.remove(name);
                                    } else {
                                        SAMRecord lookupMatch = lookup.get(name);
                                        if (lookupMatch == null) {
                                            lookup.put(name, r);
                                            continue;
                                        }
                                        if (lookupMatch.getFirstOfPairFlag() == r.getFirstOfPairFlag()) {
                                            log.error("Found multiple primary alignment records for %s", name, ". Ignoring all but first.");
                                            continue;
                                        }
                                        SAMRecord r1 = r.getFirstOfPairFlag() ? r : lookupMatch;
                                        SAMRecord r2 = r.getFirstOfPairFlag() ? lookupMatch : r;
                                        lookup.remove(name);
                                        readNames.remove(name);
                                        fq1.write(samToFastq(r1));
                                        fq2.write(samToFastq(r2));
                                    }
                                }
                            }
                        }
                    }
                    if (!lookup.isEmpty()) {
                        log.error("Missing paired primary alignment for ", readNames.size(), " reads. Writing to unpaired fastq.");
                        for (SAMRecord r : lookup.values()) {
                            fq.write(samToFastq(r));
                            readNames.remove(r.getReadName());
                        }
                    }
                }
            }
        }
        if (!readNames.isEmpty()) {
            log.error("Missing SAM records for ", readNames.size(), " reads.");
        }
    }

    public static void main(String[] argv) {
        System.exit(new ExtractFromTaxonomy().instanceMain(argv));
    }
}
