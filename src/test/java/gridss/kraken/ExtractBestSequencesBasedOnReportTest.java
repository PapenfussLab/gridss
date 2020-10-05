package gridss.kraken;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class ExtractBestSequencesBasedOnReportTest extends IntermediateFilesTest {
    @Test
    public void regression_should_extract_from_multiple_fasta() throws IOException {
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        cmd.INPUT = new File("src/test/resources/kraken/multiple_dictionaries/all_subspecies.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/multiple_dictionaries/kraken.fa"),
                new File("src/test/resources/kraken/multiple_dictionaries/virushostdb.fa")
        );
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/multiple_dictionaries/matching.nodes.dmp");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.doWork();

        Assert.assertTrue(cmd.OUTPUT.exists());
        Assert.assertTrue(cmd.REPORT_OUTPUT.exists());
        Assert.assertTrue(cmd.SUMMARY_REPORT_OUTPUT.exists());

        List<String> out = Files.readAllLines(cmd.OUTPUT.toPath());
        List<String> report = Files.readAllLines(cmd.REPORT_OUTPUT.toPath());
        List<String> summary = Files.readAllLines(cmd.SUMMARY_REPORT_OUTPUT.toPath());
        Assert.assertEquals(1, summary.size());
        Assert.assertTrue(summary.get(0).contains("papillomavirus type 45"));
        Assert.assertTrue(out.get(0).startsWith(">kraken:taxid|10593|X74479"));

    }
    @Test
    public void should_not_extract_all_strains() throws IOException {
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        cmd.INPUT = new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.kraken2.report.all.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.viral.fa")
        );
        //grep -f <(cut -f 5 < 26T.virusbreakend.vcf.kraken2.report.all.txt | sed 's/^/^/' | sed 's/$/    /') nodes.dmp > relevant_nodes.dmp
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/multistrain/relevant_nodes.dmp");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.doWork();

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
}