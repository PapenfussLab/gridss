package gridss.kraken;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class ExtractBestViralReferenceTest extends IntermediateFilesTest {
    private File fasta;
    private File report;
    private File summaryReport;
    private File intermediateTsv;
    private File finalTsv;
    private File viralReads;
    @Before
    @Override
    public void setup() throws IOException {
        super.setup();
        fasta = new File(output.toString() + ".fa");
        report = new File(output.toString() + ".report.txt");
        summaryReport = new File(output.toString() + ".summary.report.txt");
        intermediateTsv = new File(output.toString() + ".summary.tsv");
        finalTsv = new File(output.toString() + ".final.summary.tsv");
        viralReads = new File(output.toString() + ".viral.fa");
    }
    private ExtractBestViralReference setup(
            File krakenReport,
            List<File> ref,
            File nodesdmp,
            File seqidmap,
            List<Integer> taxids,
            int minSupportingReads,
            boolean favourRefSeq) {
        IdentifyViralTaxa cmd1 = new IdentifyViralTaxa();
        ExtractBestViralReference cmd2 = new ExtractBestViralReference();
        cmd1.INPUT_KRAKEN2_REPORT = krakenReport;
        cmd1.KRAKEN_REFERENCES = ref;
        cmd1.NCBI_NODES_DMP = nodesdmp;
        cmd1.SEQID2TAXID_MAP = seqidmap;
        cmd1.OUTPUT = intermediateTsv;
        cmd1.SUMMARY_REPORT_OUTPUT = summaryReport;
        cmd1.REPORT_OUTPUT = report;
        cmd1.MIN_SUPPORTING_READS = minSupportingReads;
        if (taxids != null) cmd1.TAXONOMY_IDS = taxids;
        cmd2.KRAKEN_REFERENCES = ref;
        cmd2.NCBI_NODES_DMP = nodesdmp;
        cmd2.SEQID2TAXID_MAP = seqidmap;
        cmd2.INPUT_SUMMARY = intermediateTsv;
        cmd2.OUTPUT = fasta;
        cmd2.OUTPUT_SUMMARY = finalTsv;
        cmd2.INPUT_VIRAL_READS = viralReads;
        cmd2.FAVOUR_EARLY_KRAKEN_REFERENCES = favourRefSeq;
        cmd1.doWork();
        cmd2.doWork();
        return cmd2;
    }
    @Test
    public void should_extract_from_multiple_fasta() throws IOException {
        String datadir = "src/test/resources/kraken/multiple_dictionaries/";
        viralReads.createNewFile();
        Files.write(viralReads.toPath(), "@read1\nAATACTTTTAACAATTATACTACATAAAAAAGGGTGTAACCGAAAACG\n+\nAATACTTTTAACAATTATACTACATAAAAAAGGGTGTAACCGAAAACG\n".getBytes(StandardCharsets.UTF_8));
        for (boolean prefer_refseq : new boolean[] { true, false}) {
            ExtractBestViralReference cmd = setup(
                    new File(datadir + "all_subspecies.txt"),
                    ImmutableList.of(
                            new File(datadir + "kraken.fa"),
                            new File(datadir + "virushostdb.fa")
                    ),
                    new File(datadir + "matching.nodes.dmp"),
                    new File(datadir + "seqid2taxid.map"),
                    ImmutableList.of(10239),
                    50,
                    prefer_refseq);

            List<String> out = Files.readAllLines(cmd.OUTPUT.toPath());
            Assert.assertTrue(out.get(0).startsWith(">"));
            List<String> output = Files.readAllLines(cmd.OUTPUT.toPath());
            List<String> summary = Files.readAllLines(cmd.OUTPUT_SUMMARY.toPath());
            cmd.OUTPUT.delete();
            cmd.OUTPUT_SUMMARY.delete();
            String expectedSeq = prefer_refseq ? "kraken:taxid|10593|test_in_refseq" : "kraken:taxid|10593|X74479";
            Assert.assertTrue(out.get(0).startsWith(">" + expectedSeq));
            Assert.assertTrue(summary.get(1).contains(expectedSeq));
        }
    }
    @Test
    public void should_not_extract_all_strains() throws IOException {
        String datadir = "src/test/resources/kraken/multistrain/";
        viralReads.createNewFile();
        ExtractBestViralReference cmd = setup(
                new File(datadir + "26T.virusbreakend.vcf.kraken2.report.all.txt"),
                ImmutableList.of(new File(datadir + "26T.virusbreakend.vcf.viral.fa")),
                new File(datadir + "relevant_nodes.dmp"),
                new File(datadir + "seqid2taxid.map"),
                ImmutableList.of(10239),
                1,
                false);

        Assert.assertTrue(cmd.OUTPUT.exists());
        ReferenceSequenceFile ref = new FastaSequenceFile(cmd.OUTPUT, true);
        List<ReferenceSequence> refList = new ArrayList<>();
        ReferenceSequence rs = ref.nextSequence();
        while (rs != null) {
            refList.add(rs);
            rs = ref.nextSequence();
        }
        Assert.assertEquals(1, refList.size());
    }
    @Test
    public void should_output_one_reference_per_genus() throws IOException {
        String datadir = "src/test/resources/kraken/genus/herpes_should_not_mask_HBV/";
        viralReads.createNewFile();
        ExtractBestViralReference cmd = setup(
                new File(datadir + "herpes_HBV.kraken2.report.all.txt"),
                ImmutableList.of(new File(datadir + "relevant_ref.fa")),
                new File(datadir + "relevant_nodes.dmp"),
                new File(datadir + "seqid2taxid.map"),
                ImmutableList.of(10239),
                50,
                false);
        Assert.assertTrue(cmd.OUTPUT.exists());
        ReferenceSequenceFile ref = new FastaSequenceFile(cmd.OUTPUT, true);
        List<ReferenceSequence> refList = new ArrayList<>();
        ReferenceSequence rs = ref.nextSequence();
        while (rs != null) {
            refList.add(rs);
            rs = ref.nextSequence();
        }
        Assert.assertEquals(2, refList.size());
        List<String> outputSummary = Files.readAllLines(cmd.OUTPUT_SUMMARY.toPath());
        Assert.assertEquals(2+1, outputSummary.size());
    }
}