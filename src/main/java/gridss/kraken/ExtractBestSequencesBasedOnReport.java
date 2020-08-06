package gridss.kraken;

import au.edu.wehi.idsv.kraken.KrakenReportLine;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import com.google.common.collect.Lists;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Processes a Kraken2 report and extracts the sequences with the most hits",
        oneLineSummary = "Processes a Kraken2 report and extracts the sequences with the most hits.",
        programGroup=gridss.cmdline.programgroups.DataConversion.class
)
public class ExtractBestSequencesBasedOnReport extends CommandLineProgram {
    private static final int NCBI_VIRUS_TAXID = 10239;
    private static final Log log = Log.getInstance(ExtractBestSequencesBasedOnReport.class);
    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Kraken2 report file.")
    public File INPUT;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output fastq file")
    public File OUTPUT;
    @Argument(doc="NCBI Taxonomy IDs to extract. All taxonomic entries under these IDs are also extracted. Defaults to all viruses.")
    public List<Integer> TAXONOMY_IDS = Lists.newArrayList(NCBI_VIRUS_TAXID);
    @Argument(doc="NCBI taxonomy nodes.dmp. Download and extract from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
    public File NCBI_NODES_DMP;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsReadable(NCBI_NODES_DMP);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            IndexedFastaSequenceFile ref = new IndexedFastaSequenceFile(REFERENCE_SEQUENCE);
            log.info("Parsing Kraken2 report from ", INPUT);
            List<KrakenReportLine> report = Files.lines(INPUT.toPath())
                    .map(s -> new KrakenReportLine(s))
                    .collect(Collectors.toList());
            log.info("Loading NCBI taxonomy from ", NCBI_NODES_DMP);
            boolean[] taxIdLookup = TaxonomyHelper.createInclusionLookup(TAXONOMY_IDS, TaxonomyHelper.parseMinimal(NCBI_NODES_DMP));
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    public static void main(String[] argv) {
        System.exit(new ExtractBestSequencesBasedOnReport().instanceMain(argv));
    }
}
