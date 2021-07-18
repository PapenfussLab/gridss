package gridss.kraken;

import au.edu.wehi.idsv.debruijn.ContigKmerCounter;
import au.edu.wehi.idsv.kraken.KrakenReportLine;
import au.edu.wehi.idsv.kraken.SeqIdToTaxIdMap;
import au.edu.wehi.idsv.ncbi.MinimalTaxonomyNode;
import au.edu.wehi.idsv.ncbi.TaxonomyHelper;
import au.edu.wehi.idsv.ncbi.TaxonomyLevel;
import au.edu.wehi.idsv.ncbi.TaxonomyNode;
import com.google.common.collect.Streams;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMSequenceDictionary;
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
import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang.StringUtils;
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
import java.util.stream.Stream;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

@CommandLineProgramProperties(
        summary = "Processes a Kraken2 report and extracts the sequences with the most hits",
        oneLineSummary = "Processes a Kraken2 report and extracts the sequences with the most hits.",
        programGroup=gridss.cmdline.programgroups.DataConversion.class
)
public class ExtractBestViralReference extends CommandLineProgram {
    private static final Log log = Log.getInstance(ExtractBestViralReference.class);
    private static final Comparator<KrakenReportLine> SORT_ORDER = KrakenReportLine.ByCountAssignedDirectly.reversed().thenComparing(KrakenReportLine.ByCountAssignedToTree.reversed());
    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="TSV from gridss.IdentifyViralTaxa")
    public File INPUT_SUMMARY;
    @Argument(doc="Viral reads. Used to determine which genome for the chosen taxa to return.")
    public File INPUT_VIRAL_READS;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output fasta file")
    public File OUTPUT;
    @Argument(doc="TSV annotated with extracted references")
    public File OUTPUT_SUMMARY;
    @Argument(doc="TSV containing number of matched kmers for each candidate reference sequence.")
    public File OUTPUT_MATCHING_KMERS;
    @Argument(doc="Kmer used determining best viral genome match for viral reads.")
    public int KMER = 16;
    @Argument(doc="Distance between kmers in reference lookup. Longer stride reduces memory usage. Should not be more than kmer length")
    public int STRIDE = 16;
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
    @Argument(doc="Maximum number of contigs to extract per NCBI taxonomic identifiers.", optional = true)
    public int CONTIGS_PER_TAXID = 1;
    @Argument(doc="Use the viral references that occur in earlier KRAKEN_REFERENCES files whenever possible.", optional = true)
    public boolean FAVOUR_EARLY_KRAKEN_REFERENCES = true;

    @Override
    protected String[] customCommandLineValidation() {
        if (KRAKEN_REFERENCES == null || KRAKEN_REFERENCES.size() == 0) {
            return new String[] {"KRAKEN_REFERENCES required. At minimum, this is the library/viral/library.fna file in the kraken2 database directory." };
        }
        return super.customCommandLineValidation();
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(SEQID2TAXID_MAP);
        IOUtil.assertFileIsReadable(NCBI_NODES_DMP);
        IOUtil.assertFileIsReadable(INPUT_SUMMARY);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (OUTPUT_SUMMARY != null) IOUtil.assertFileIsWritable(OUTPUT_SUMMARY);
        try {
            List<IndexedFastaSequenceFile> ref = new ArrayList<>(KRAKEN_REFERENCES.size());
            for (File f : KRAKEN_REFERENCES) {
                IOUtil.assertFileIsReadable(f);
                ReferenceCommandLineProgram.ensureSequenceDictionary(f);
                ref.add(new IndexedFastaSequenceFile(f));
            }
            log.info("Loading seqid2taxid.map from ", SEQID2TAXID_MAP);
            Map<String, Integer> seq2taxLookup = SeqIdToTaxIdMap.createLookup(SEQID2TAXID_MAP);
            log.info("Loading NCBI taxonomy from ", NCBI_NODES_DMP);
            Map<Integer, TaxonomyNode> taxa = TaxonomyHelper.parseFull(NCBI_NODES_DMP);
            log.info("Parsing ", INPUT_SUMMARY);
            List<List<String>> summaryLines = Files.readAllLines(INPUT_SUMMARY.toPath()).stream()
                    .map(line -> Arrays.asList(line.split("\t")))
                    .collect(Collectors.toList());
            Set<Integer> taxaOfInterest = summaryLines.stream()
                    .skip(1) // ignore header
                    .map(line -> Integer.parseInt(line.get(6)))
                    .collect(Collectors.toSet());
            boolean[] taxidInTreeOfInterest = TaxonomyHelper.createInclusionLookup(taxaOfInterest, taxa);

            Stream<ReferenceSequence> candidateContigs = ref.stream()
                .flatMap(r -> r.getSequenceDictionary()
                    .getSequences()
                    .stream()
                    .filter(s -> taxidInTreeOfInterest[seq2taxLookup.get(s.getSequenceName())])
                    .map(s -> r.getSequence(s.getSequenceName())));

            ContigKmerCounter ckc = new ContigKmerCounter(candidateContigs, KMER, STRIDE);
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
                .collect(groupingBy(p -> getParentTaxaOfInterest(taxaOfInterest, taxa, seq2taxLookup.get(p.getKey()))));
            if (OUTPUT_MATCHING_KMERS != null) {
                log.info("Writing matching kmer counts to ", OUTPUT_MATCHING_KMERS);
                Files.write(OUTPUT_MATCHING_KMERS.toPath(), candidateContigCountsByTaxa.keySet()
                    .stream()
                    .flatMap(parentTaxId -> candidateContigCountsByTaxa.get(parentTaxId)
                            .stream()
                            .map(p -> String.format("%d\t%s\t%d", parentTaxId, p.getKey(), p.getValue())))
                    .collect(Collectors.toList()));
            }
            List<List<String>> output = new ArrayList<>();
            List<String> outputReferences = new ArrayList<>();
            for (int i = 0; i < summaryLines.size(); i++) {
                List<String> line = summaryLines.get(i);
                if (i == 0) {
                    List<String> outputLine = new ArrayList<>(line);
                    outputLine.add("reference");
                    //outputLine.add("reference_tax_name"); // Need names.dmp to get this since it's not necessarily in the kraken2 report
                    outputLine.add("reference_taxid");
                    outputLine.add("reference_kmer_count");
                    outputLine.add("alternate_kmer_count");
                    output.add(outputLine);
                } else {
                    int kraken2taxid = Integer.parseInt(line.get(6));
                    Map<String, Integer> krakenReferenceFileIndex = new HashMap<>();
                    for (int krOffset = 0; krOffset < ref.size(); krOffset++) {
                        IndexedFastaSequenceFile r = ref.get(krOffset);
                        for (SAMSequenceRecord sr : r.getSequenceDictionary().getSequences()) {

                        }
                    }
                    List<Pair<String, Long>> candidatesForTaxa = candidateContigCountsByTaxa.get(kraken2taxid).stream()
                            .sorted(
                                    Comparator.comparingInt((Pair<String, Long> p) -> FAVOUR_EARLY_KRAKEN_REFERENCES ? offsetInList(ref, p.getKey()) : 0)
                                    .thenComparing(Comparator.comparingLong((Pair<String, Long> p) -> p.getValue()).reversed()))
                            .collect(Collectors.toList());
                    List<Pair<String, Long>> outputRef = candidatesForTaxa.stream().limit(CONTIGS_PER_TAXID).collect(Collectors.toList());
                    List<Pair<String, Long>> missingRef = candidatesForTaxa.stream().skip(CONTIGS_PER_TAXID).collect(Collectors.toList());
                    log.debug("Processing taxid " + kraken2taxid);
                    for (Pair<String, Long> pair : candidatesForTaxa) {
                        log.debug(pair.getKey() + "\t" + pair.getValue());
                    }
                    long nextBestKmerCount = missingRef.stream().mapToLong(p -> p.getValue()).findFirst().orElse(0);
                    for (Pair<String, Long> pair : outputRef) {
                        String refName = pair.getKey();
                        Long kmersAssigned = pair.getValue();
                        int referenceTaxId = seq2taxLookup.get(refName);
                        log.info("Using " + pair.getKey() + " as viral reference for " + kraken2taxid + " (" + referenceTaxId + ")");
                        List<String> outputLine = new ArrayList<>(line);
                        outputLine.add(refName);
                        //outputLine.add(taxa.get(referenceTaxId).name); // name isn't in nodes.dmp
                        outputLine.add(Integer.toString(referenceTaxId));
                        outputLine.add(kmersAssigned.toString());
                        outputLine.add(Long.toString(nextBestKmerCount));
                        output.add(outputLine);
                        outputReferences.add(refName);
                    }
                }
            }
            log.info("Writing " + outputReferences.size() + " viral contigs to ", OUTPUT);
            if (outputReferences.size() > 0) {
                try (FastaReferenceWriter writer = new FastaReferenceWriterBuilder()
                        .setMakeDictOutput(true)
                        .setMakeFaiOutput(true)
                        .setFastaFile(OUTPUT.toPath())
                        .build()) {
                    for (String refContig : outputReferences) {
                        for (IndexedFastaSequenceFile fa : ref) {
                            if (fa.getSequenceDictionary().getSequence(refContig) != null) {
                                ReferenceSequence seq = fa.getSequence(refContig);
                                writer.addSequence(cleanSequence(seq));
                                break;
                            }
                        }
                    }
                }
            } else {
                // workaround for https://github.com/samtools/htsjdk/issues/1498
                Files.write(OUTPUT.toPath(), new byte[0]);
                log.warn("No sequences written to ", OUTPUT);
            }
            List<String> finalOutput = output.stream()
                    .map(line -> StringUtils.join(line, "\t"))
                    .collect(Collectors.toList());
            Files.write(OUTPUT_SUMMARY.toPath(), finalOutput);
        } catch (IOException e) {
            log.error(e);
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    private int getParentTaxaOfInterest(Set<Integer> taxaOfInterest, Map<Integer, ? extends MinimalTaxonomyNode> taxa, int taxId) {
        while (!taxaOfInterest.contains(taxId)) {
            int parentTaxId = taxa.get(taxId).parentTaxId;
            if (parentTaxId <= 1) return taxId;
            taxId = parentTaxId;
        }
        return taxId;
    }

    private static int offsetInList(List<IndexedFastaSequenceFile> ref, String contigName) {
        for (int i = 0; i < ref.size(); i++) {
            if (ref.get(i).getSequenceDictionary().getSequence(contigName) != null) {
                return i;
            }
        }
        return ref.size();
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

    public static void main(String[] argv) {
        System.exit(new ExtractBestViralReference().instanceMain(argv));
    }
}
