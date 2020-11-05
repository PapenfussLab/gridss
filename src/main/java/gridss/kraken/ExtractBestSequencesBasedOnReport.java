package gridss.kraken;

import au.edu.wehi.idsv.kraken.KrakenReportLine;
import au.edu.wehi.idsv.ncbi.MinimalTaxonomyNode;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import au.edu.wehi.idsv.ncbi.TaxonomyLevel;
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
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.groupingBy;

@CommandLineProgramProperties(
        summary = "Processes a Kraken2 report and extracts the sequences with the most hits",
        oneLineSummary = "Processes a Kraken2 report and extracts the sequences with the most hits.",
        programGroup=gridss.cmdline.programgroups.DataConversion.class
)
public class ExtractBestSequencesBasedOnReport extends CommandLineProgram {
    private static final int NCBI_VIRUS_TAXID = 10239;
    private static final Log log = Log.getInstance(ExtractBestSequencesBasedOnReport.class);
    private static final Comparator<KrakenReportLine> SORT_ORDER = KrakenReportLine.ByCountAssignedDirectly.reversed().thenComparing(KrakenReportLine.ByCountAssignedToTree.reversed());
    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Kraken2 report file.")
    public File INPUT;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output fasta file")
    public File OUTPUT;
    @Argument(shortName="RO", doc="Kraken2 report filtered to only taxa of interest.", optional = true)
    public File REPORT_OUTPUT;
    @Argument(shortName="SRO", doc="Kraken2 report filtered to only taxa included in the output fasta file.", optional = true)
    public File SUMMARY_REPORT_OUTPUT;
    @Argument(shortName="SO", doc="Summary csv", optional = true)
    public File SUMMARY_OUTPUT;
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
    @Argument(doc="Maximum number of NCBI taxonomic identifiers to extract sequences for.", optional = true)
    public Integer TAXA_TO_RETURN = null;
    @Argument(doc="Maximum number of contigs to extract per NCBI taxonomic identifiers.", optional = true)
    public int CONTIGS_PER_TAXID = 1;
    @Argument(doc="Minimum number of supporting reads", optional = true)
    public int MIN_SUPPORTING_READS = 1;
    @Argument(doc="Taxonomic level for which only one sequence will be output. Useful to prevent multiple strains of the same/similar viruses being output.", optional = true)
    public TaxonomyLevel TAXONOMIC_DEDUPLICATION_LEVEL = TaxonomyLevel.Genus;
    public int TAXA_PER_DEDUPLICATION_LEVEL = 1;

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
        if (SUMMARY_REPORT_OUTPUT != null) IOUtil.assertFileIsWritable(SUMMARY_REPORT_OUTPUT);
        if (OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);
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
            List<KrakenReportLine> fullReport = Files.lines(INPUT.toPath())
                    .map(s -> new KrakenReportLine(s))
                    .collect(Collectors.toList());
            Int2IntMap taxaGroupLookup = createTaxaGroupLookup(taxa, fullReport, TAXONOMIC_DEDUPLICATION_LEVEL);
            List<KrakenReportLine> filteredReport = fullReport.stream()
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
            Map<Integer, List<KrakenReportLine>> interestingGroups = filteredReport.stream()
                    .filter(s -> taxaWithSequence[s.taxonomyId])
                    .filter(s -> s.countAssignedToTree >= MIN_SUPPORTING_READS)
                    .collect(groupingBy(krl -> taxaGroupLookup.getOrDefault(krl.taxonomyId, krl.taxonomyId)));
            List<KrakenReportLine> interestingNodes = interestingGroups.values().stream()
                    .flatMap(list -> list.stream().sorted(SORT_ORDER).limit(TAXA_PER_DEDUPLICATION_LEVEL))
                    .sorted(SORT_ORDER)
                    .limit(TAXA_TO_RETURN == null ? Integer.MAX_VALUE : TAXA_TO_RETURN)
                    .collect(Collectors.toList());
            if (SUMMARY_REPORT_OUTPUT != null) {
                log.info("Writing summary kraken report to ", SUMMARY_REPORT_OUTPUT);
                Files.write(SUMMARY_REPORT_OUTPUT.toPath(), interestingNodes.stream().map(krl -> krl.line).collect(Collectors.toList()));
            }
            List<List<Integer>> sequenceIndexesToExport = new ArrayList<>();
            for (IndexedFastaSequenceFile fa : ref) {
                sequenceIndexesToExport.add(new ArrayList<>());
            }
            Map<KrakenReportLine, List<String>> exportedRefs = new HashMap<>();
            for (KrakenReportLine taxaToExport : interestingNodes) {
                List<String> refsForThisTaxa = new ArrayList<>();
                int taxid = taxaToExport.taxonomyId;
                int contigsFound = 0;
                for (int i = 0; i < ref.size(); i++) {
                    IndexedFastaSequenceFile fa = ref.get(i);
                    for (SAMSequenceRecord ssr : fa.getSequenceDictionary().getSequences()) {
                        int seqTaxId = extractTaxIdFromKrakenSequence(ssr);
                        if (seqTaxId == taxid && contigsFound < CONTIGS_PER_TAXID) {
                            contigsFound++;
                            sequenceIndexesToExport.get(i).add(ssr.getSequenceIndex());
                        }
                    }
                }
                if (contigsFound > 0) {
                    exportedRefs.put(taxaToExport, refsForThisTaxa);
                }
            }
            if (SUMMARY_OUTPUT != null) {
                log.info("Writing summary csv report to ", SUMMARY_OUTPUT);
                List<String> summary = new ArrayList<>();
                summary.add(createSummaryHeader());
                for (Map.Entry<KrakenReportLine, List<String>> in : exportedRefs.entrySet()) {
                    summary.add(createSummaryLine(fullReport, taxa, in.getKey(), in.getValue()));
                }
                Files.write(SUMMARY_OUTPUT.toPath(), summary);
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

    private String createSummaryHeader() {
        return "taxid_genus\tname_genus\treads_genus\ttaxid_species\treads_species\tname_species\ttaxid\tname\treads\treference";
    }

    private String createSummaryLine(List<KrakenReportLine> fullReport, Map<Integer, MinimalTaxonomyNode> taxa, KrakenReportLine line, List<String> ref) {
        Map<Integer, KrakenReportLine> lookup = fullReport.stream().collect(Collectors.toMap(x -> x.taxonomyId, x -> x));
        KrakenReportLine genus = line;
        KrakenReportLine species = line;
        KrakenReportLine current = line;
        while (current != null && current.taxonomyId > 1) {
            switch (current.rank) {
                case "S":
                    species = line;
                    break;
                case "G":
                    genus = line;
                    break;
            }
            int parent_taxid = taxa.get(current.taxonomyId).parentTaxId;
            if (parent_taxid <= 1) break;
            current = lookup.get(parent_taxid);
        }
        return String.format("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s",
                genus.taxonomyId, genus.scientificName.trim(),genus.countAssignedToTree,
                genus.taxonomyId, species.scientificName.trim(),species.countAssignedToTree,
                line.taxonomyId, line.scientificName.trim(),line.countAssignedToTree,
                ref.stream().collect(Collectors.joining(",")));
    }

    /**
     * Uses the kraken report to create a lookup at the given taxonomic level.
     * @param fullReport
     * @param level
     * @return
     */
    public static Int2IntMap createTaxaGroupLookup(Map<Integer, MinimalTaxonomyNode> taxa, List<KrakenReportLine> fullReport, TaxonomyLevel level) {
        Map<Integer, KrakenReportLine> reportLookup = fullReport.stream().collect(Collectors.toMap(krl -> krl.taxonomyId, krl -> krl));
        Int2IntMap taxaGroupLookup = new Int2IntOpenHashMap();
        for (KrakenReportLine krl : fullReport) {
            int groupTaxa = krl.taxonomyId;
            if (level != null) {
                KrakenReportLine parent = krl;
                while (parent != null) {
                    if (level.krakenAbbreviation().equals(parent.rank)) {
                        groupTaxa = parent.taxonomyId;
                    }
                    if (taxa.get(parent.taxonomyId) == null) break;
                    if (parent.taxonomyId == taxa.get(parent.taxonomyId).parentTaxId) break;
                    parent = reportLookup.get(taxa.get(parent.taxonomyId).parentTaxId);
                }
            }
            taxaGroupLookup.put(krl.taxonomyId, groupTaxa);
        }
        return taxaGroupLookup;
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
