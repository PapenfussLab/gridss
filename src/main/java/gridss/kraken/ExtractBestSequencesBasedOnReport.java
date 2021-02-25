package gridss.kraken;

import au.edu.wehi.idsv.debruijn.ContigKmerCounter;
import au.edu.wehi.idsv.kraken.KrakenReportLine;
import au.edu.wehi.idsv.ncbi.MinimalTaxonomyNode;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import au.edu.wehi.idsv.ncbi.TaxonomyLevel;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Streams;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
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
import org.apache.commons.lang3.tuple.Pair;
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
    public File INPUT_KRAKEN2_REPORT;
    @Argument(doc="Viral reads. Used to determine which genome for the chosen taxa to return.")
    public File INPUT_VIRAL_READS;
    @Argument(doc="Kmer used determining best viral genome match for viral reads.")
    public int VIRAL_READ_KMER = 31;
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
    @Argument(doc="Kraken2 seqid2taxid.map mapping file")
    public File SEQID2TAXID_MAP;
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
        IOUtil.assertFileIsReadable(INPUT_KRAKEN2_REPORT);
        IOUtil.assertFileIsReadable(NCBI_NODES_DMP);
        IOUtil.assertFileIsReadable(SEQID2TAXID_MAP);
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
            log.info("Loading seqid2taxid.map from ", SEQID2TAXID_MAP);
            Map<String, Integer> seq2taxLookup = createSeqId2TaxIdMap(SEQID2TAXID_MAP);
            log.info("Loading NCBI taxonomy from ", NCBI_NODES_DMP);
            Map<Integer, MinimalTaxonomyNode> taxa = TaxonomyHelper.parseMinimal(NCBI_NODES_DMP);
            boolean[] taxIdLookup = TaxonomyHelper.createInclusionLookup(TAXONOMY_IDS, taxa);
            boolean[] relevantTaxIdAndAncestors = TaxonomyHelper.addAncestors(taxIdLookup, taxa);
            log.info("Parsing Kraken2 report from ", INPUT_KRAKEN2_REPORT);
            List<KrakenReportLine> fullReport = Files.lines(INPUT_KRAKEN2_REPORT.toPath())
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
                    .mapToInt(s -> seq2taxLookup.get(s.getSequenceName()))
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
            Set<Integer> taxaOfInterest = interestingNodes.stream()
                    .map(krl -> krl.taxonomyId)
                    .collect(Collectors.toSet());
            List<ReferenceSequence> candidateContigs = ref.stream()
                .flatMap(r -> r.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(s -> taxaOfInterest.contains(seq2taxLookup.get(s.getSequenceName())))
                    .map(s -> r.getSequence(s.getSequenceName())))
                .collect(Collectors.toList());
            ContigKmerCounter ckc = new ContigKmerCounter(candidateContigs, VIRAL_READ_KMER);
            if (INPUT_VIRAL_READS != null) {
                log.info("Identifying best viral reference genomes from ", INPUT_VIRAL_READS);
                FastqReader fqr = new FastqReader(INPUT_VIRAL_READS);
                while (fqr.hasNext()) {
                    FastqRecord fq = fqr.next();
                    ckc.count(fq.getReadBases());
                }
            }
            Map<Integer, List<Pair<String, Long>>> candidateContigCountsByTaxa = Streams.zip(
                    ckc.getContigs().stream(),
                    ckc.getKmerCounts().stream(),
                    Pair::of)
                .collect(groupingBy(p -> seq2taxLookup.get(p.getKey())));
            Set<String> contigsInViralReference = new HashSet<>();
            for (KrakenReportLine krl : interestingNodes) {
                List<Pair<String, Long>> candidatesForTaxa = candidateContigCountsByTaxa.get(krl.taxonomyId).stream()
                        .sorted(Comparator.comparing((Pair<String, Long> p) -> p.getValue()).reversed())
                        .collect(Collectors.toList());
                for (Pair<String, Long> p : candidatesForTaxa) {
                    log.debug("taxid\t" + krl.taxonomyId + "\t" + p.getKey() + "\t" + p.getValue());
                }
                contigsInViralReference.addAll(
                    candidatesForTaxa.stream()
                    .limit(CONTIGS_PER_TAXID)
                    .map(p -> p.getKey())
                    .collect(Collectors.toSet()));
            }
            List<String> summary = new ArrayList<>();
            summary.add(createSummaryHeader());
            log.info("Writing " + contigsInViralReference.size() + " viral contigs to ", OUTPUT);
            if (contigsInViralReference.size() > 0) {
                try (FastaReferenceWriter writer = new FastaReferenceWriterBuilder()
                        .setMakeDictOutput(true)
                        .setMakeFaiOutput(true)
                        .setFastaFile(OUTPUT.toPath())
                        .build()) {
                    for (IndexedFastaSequenceFile fa : ref) {
                        for (SAMSequenceRecord ssr : fa.getSequenceDictionary().getSequences()) {
                            if (contigsInViralReference.contains(ssr.getSequenceName())) {
                                // only export the first contig if we find the same contig name in multiple files
                                contigsInViralReference.remove(ssr.getSequenceName());
                                writer.addSequence(cleanSequence(fa.getSequence(ssr.getSequenceName())));
                                int taxid = seq2taxLookup.get(ssr.getSequenceName());
                                KrakenReportLine krl = interestingNodes.stream().filter(x -> x.taxonomyId == taxid).findFirst().get();
                                summary.add(createSummaryLine(fullReport, taxa, krl, ssr.getSequenceName()));
                            }
                        }
                    }
                }
            } else {
                // workaround for https://github.com/samtools/htsjdk/issues/1498
                Files.write(OUTPUT.toPath(), new byte[0]);
                log.warn("No sequences written to ", OUTPUT);
            }
            if (SUMMARY_OUTPUT != null) {
                log.info("Writing summary csv report to ", SUMMARY_OUTPUT);
                Files.write(SUMMARY_OUTPUT.toPath(), summary);
            }
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    private String createSummaryHeader() {
        return "taxid_genus\tname_genus\treads_genus\ttaxid_species\tname_species\treads_species\ttaxid\tname\treads\treference";
    }

    private String createSummaryLine(List<KrakenReportLine> fullReport, Map<Integer, MinimalTaxonomyNode> taxa, KrakenReportLine line, String ref) {
        Map<Integer, KrakenReportLine> lookup = fullReport.stream().collect(Collectors.toMap(x -> x.taxonomyId, x -> x));
        KrakenReportLine genus = line;
        KrakenReportLine species = line;
        KrakenReportLine current = line;
        while (current != null && current.taxonomyId > 1) {
            switch (current.rank) {
                case "S":
                    species = current;
                    break;
                case "G":
                    genus = current;
                    break;
            }
            int parent_taxid = taxa.get(current.taxonomyId).parentTaxId;
            if (parent_taxid <= 1) break;
            current = lookup.get(parent_taxid);
        }
        return String.format("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s",
                genus.taxonomyId, genus.scientificName.trim(),genus.countAssignedToTree,
                species.taxonomyId, species.scientificName.trim(),species.countAssignedToTree,
                line.taxonomyId, line.scientificName.trim(),line.countAssignedToTree,
                ref);
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
    private static Map<String, Integer> createSeqId2TaxIdMap(File seqid2taxid_map) {
        try {
            List<String> lines = Files.readAllLines(seqid2taxid_map.toPath());
            return lines.stream()
                .map(s -> s.split("\t"))
                .collect(Collectors.toMap(s -> s[0], s -> Integer.parseInt(s[1])));
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    public static void main(String[] argv) {
        System.exit(new ExtractBestSequencesBasedOnReport().instanceMain(argv));
    }
}
