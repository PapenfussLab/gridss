package gridss.kraken;

import au.edu.wehi.idsv.kraken.KrakenReportLine;
import au.edu.wehi.idsv.ncbi.MinimalTaxonomyNode;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import com.google.common.collect.Lists;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
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
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output fasta file")
    public File OUTPUT;
    @Argument(shortName="RO", doc="Kraken2 report filtered to only taxa of interest.", optional = true)
    public File REPORT_OUTPUT;
    @Argument(shortName="SRO", doc="Kraken2 report filtered to only taxa included in the output fasta file.", optional = true)
    public File SUMMARY_REPORT_OUTPUT;
    @Argument(doc="NCBI Taxonomy IDs to extract. All taxonomic entries under these IDs are also extracted. Defaults to all viruses.")
    public List<Integer> TAXONOMY_IDS = Lists.newArrayList(NCBI_VIRUS_TAXID);
    @Argument(doc="NCBI taxonomy nodes.dmp. Download and extract from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip")
    public File NCBI_NODES_DMP;
    @Argument(doc="Kraken2 library.fna files." +
            " Downloaded by kraken2-build." +
            " Must be indexed." +
            " Do not run kraken2-build --clean as these files will be removed." +
            " Files are checked in order and all the contigs for the given taxid from the first matching file are extracted.", optional = true)
    public List<File> KRAKEN_REFERENCES;
    @Argument(doc="Maximum number of taxid to extract sequences for.", optional = true)
    public int TAXID_TO_RETURN = 1;
    @Argument(doc="Minimum number of supporting reads", optional = true)
    public int MIN_SUPPORTING_READS = 1;

    @Override
    protected String[] customCommandLineValidation() {
        if (KRAKEN_REFERENCES == null || KRAKEN_REFERENCES.size() == 0) {
            return new String[] {"KRAKEN_REFERENCES required. This file is located under library/viral/library.fna in the kraken2 database directory." };
        }
        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(NCBI_NODES_DMP);
        IOUtil.assertFileIsWritable(OUTPUT);
        try {
            List<IndexedFastaSequenceFile> ref = new ArrayList<>(KRAKEN_REFERENCES.size());
            for (File f : KRAKEN_REFERENCES) {
                IOUtil.assertFileIsReadable(f);
                ReferenceCommandLineProgram.ensureSequenceDictionary(f);
                ref.add(new IndexedFastaSequenceFile(f));
            }
            log.info("Loading NCBI taxonomy from ", NCBI_NODES_DMP);
            Map<Integer, MinimalTaxonomyNode> taxa = TaxonomyHelper.parseMinimal(NCBI_NODES_DMP);
            boolean[] taxIdLookup = TaxonomyHelper.createInclusionLookup(TAXONOMY_IDS, taxa);
            boolean[] relevantTaxIdAndAncestors = TaxonomyHelper.addAncestors(taxIdLookup, taxa);
            log.info("Parsing Kraken2 report from ", INPUT);
            List<KrakenReportLine> filteredReport = Files.lines(INPUT.toPath())
                    .map(s -> new KrakenReportLine(s))
                    .filter(s -> relevantTaxIdAndAncestors[s.taxonomyId])
                    .collect(Collectors.toList());
            if (REPORT_OUTPUT != null) {
                log.info("Writing abridged report to ", REPORT_OUTPUT);
                Files.write(REPORT_OUTPUT.toPath(), filteredReport.stream().map(krl -> krl.line).collect(Collectors.toList()));
            }
            boolean[] taxaWithSequence = new boolean[taxIdLookup.length];
            ref.stream()
                    .flatMap(r -> r.getSequenceDictionary().getSequences().stream())
                    .mapToInt(s -> extractTaxIdFromKrakenSequence(s))
                    .forEach(x -> taxaWithSequence[x] = true);
            List<KrakenReportLine> interestingNodes = filteredReport.stream()
                    .filter(s -> taxaWithSequence[s.taxonomyId])
                    .filter(s -> s.countAssignedToTree >= MIN_SUPPORTING_READS)
                    .sorted(KrakenReportLine.ByCountAssignedDirectly.reversed().thenComparing(KrakenReportLine.ByCountAssignedToTree.reversed()))
                    .limit(TAXID_TO_RETURN)
                    .collect(Collectors.toList());
            if (SUMMARY_REPORT_OUTPUT != null) {
                log.info("Writing summary report to ", SUMMARY_REPORT_OUTPUT);
                Files.write(SUMMARY_REPORT_OUTPUT.toPath(), interestingNodes.stream().map(krl -> krl.line).collect(Collectors.toList()));
            }
            List<List<Integer>> sequenceIndexesToExport = new ArrayList<>();
            for (IndexedFastaSequenceFile fa : ref) {
                sequenceIndexesToExport.add(new ArrayList<>());
            }
            for (KrakenReportLine taxaToExport : interestingNodes) {
                int taxid = taxaToExport.taxonomyId;
                boolean found = false;
                for (int i = 0; i < ref.size() && !found; i++) {
                    IndexedFastaSequenceFile fa = ref.get(i);
                    for (SAMSequenceRecord ssr : fa.getSequenceDictionary().getSequences()) {
                        int seqTaxId = extractTaxIdFromKrakenSequence(ssr);
                        if (seqTaxId == taxid) {
                            sequenceIndexesToExport.get(i).add(ssr.getSequenceIndex());
                            found = true;
                        }
                    }
                }
            }
            if (sequenceIndexesToExport.stream().anyMatch(list -> !list.isEmpty())) {
                try (FastaReferenceWriter writer = new FastaReferenceWriterBuilder()
                        .setMakeDictOutput(true)
                        .setMakeFaiOutput(true)
                        .setFastaFile(OUTPUT.toPath())
                        .build()) {
                    for (int i = 0; i < ref.size(); i++) {
                        IndexedFastaSequenceFile fa = ref.get(i);
                        List<Integer> indexes = sequenceIndexesToExport.get(i);
                        Collections.sort(indexes);
                        for (int j : indexes) {
                            SAMSequenceRecord seq = fa.getSequenceDictionary().getSequence(j);
                            writer.addSequence(cleanSequence(fa.getSequence(seq.getSequenceName())));
                        }
                    }
                }
            } else {
                // workaround for https://github.com/samtools/htsjdk/issues/1498
                Files.write(OUTPUT.toPath(), new byte[0]);
                log.warn("No sequences written to ", OUTPUT);
            }
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    /**
     * Some sequences contain invalid/masked bases. These are converted to Ns.
     * @param rs
     * @return
     */
    private static ReferenceSequence cleanSequence(ReferenceSequence rs) {
        byte[] seq = rs.getBases();
        for (int i = 0; i < seq.length; i++) {
            if (!SequenceUtil.isIUPAC(seq[i])) {
                seq[i] = 'N';
            }
        }
        return new ReferenceSequence(rs.getName(), rs.getContigIndex(), seq);
    }
    private static int extractTaxIdFromKrakenSequence(SAMSequenceRecord ssr) {
        return Integer.parseInt(ssr.getSequenceName().split("[|]")[1]);
    }

    public static void main(String[] argv) {
        System.exit(new ExtractBestSequencesBasedOnReport().instanceMain(argv));
    }
}
