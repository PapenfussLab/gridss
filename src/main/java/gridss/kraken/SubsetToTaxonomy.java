package gridss.kraken;

import au.edu.wehi.idsv.kraken.KrakenClassification;
import au.edu.wehi.idsv.kraken.KrakenParser;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Processes Kraken2 output and subsets to only those records under the given taxonomic IDs. ",
        oneLineSummary = "Processes Kraken2 output and subsets to only those records under the given taxonomic IDs.",
        programGroup=gridss.cmdline.programgroups.DataConversion.class
)
public class SubsetToTaxonomy extends CommandLineProgram {
    private static final int NCBI_VIRUS_TAXID = 10239;
    private static final Log log = Log.getInstance(SubsetToTaxonomy.class);
    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Kraken2 output file.")
    public File INPUT;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file")
    public File OUTPUT;
    @Argument(shortName = "F", doc="Output format. Valid values are KRAKEN, and READ_NAME", optional = true)
    public OutputFormat FORMAT = OutputFormat.KRAKEN;
    @Argument(doc="NCBI Taxonomy IDs to extract. All taxonomic entries under these IDs are also extracted. Defaults to all viruses.")
    public List<Integer> TAXONOMY_IDS = Lists.newArrayList(NCBI_VIRUS_TAXID);
    @Argument(doc="NCBI taxonomy nodes.dmp. Download and extract from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
    public File NCBI_NODES_DMP;

    public enum OutputFormat {
        /**
         * Pass-through of Kraken2 output
         */
        KRAKEN,
        /**
         * Read name
         */
        READ_NAME,
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(NCBI_NODES_DMP);
        IOUtil.assertFileIsWritable(OUTPUT);
        try (KrakenParser parser = new KrakenParser(new BufferedReader(new InputStreamReader(new FileInputStream(INPUT))))) {
            log.info("Loading NCBI taxonomy from ", NCBI_NODES_DMP);
            boolean[] taxIdLookup = TaxonomyHelper.createInclusionLookup(TAXONOMY_IDS, TaxonomyHelper.parseMinimal(NCBI_NODES_DMP));
            log.info("Performing taxonomy lookup on ", INPUT);
            try (BufferedOutputStream os = new BufferedOutputStream(new FileOutputStream(OUTPUT))) {
                while (parser.hasNext()) {
                    KrakenClassification kc = parser.next();
                    if (taxIdLookup[kc.taxonomyId]) {
                        switch (FORMAT) {
                            case READ_NAME:
                                os.write(kc.sequenceId.getBytes(StandardCharsets.UTF_8));
                                os.write('\n');
                            case KRAKEN:
                            default:
                                os.write(kc.toKrakenOutput().getBytes(StandardCharsets.UTF_8));
                                os.write('\n');
                        }
                    }
                }
            }
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    public static void main(String[] argv) {
        System.exit(new SubsetToTaxonomy().instanceMain(argv));
    }
}
