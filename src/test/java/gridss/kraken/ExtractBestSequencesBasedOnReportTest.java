package gridss.kraken;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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

}