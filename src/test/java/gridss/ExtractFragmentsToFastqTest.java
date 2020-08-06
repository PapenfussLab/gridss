package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;

public class ExtractFragmentsToFastqTest extends IntermediateFilesTest {
    private List<List<FastqRecord>> go(Collection<String> readNames, SAMRecord... records) throws IOException {
        File rnFile = new File(testFolder.getRoot(), "readnames.txt");
        File fq = new File(testFolder.getRoot(), "fq.fq");
        File fq1 = new File(testFolder.getRoot(), "fq1.fq");
        File fq2 = new File(testFolder.getRoot(), "fq2.fq");
        Files.write(rnFile.toPath(), readNames);
        createInput(records);
        ExtractFragmentsToFastq cmd = new ExtractFragmentsToFastq();
        cmd.INPUT = input;
        cmd.OUTPUT_FQ = fq;
        cmd.OUTPUT_FQ1 = fq1;
        cmd.OUTPUT_FQ2 = fq2;
        cmd.READ_NAMES = rnFile;
        cmd.doWork();
        List<List<FastqRecord>> result = ImmutableList.of(
                getFastqRecords(fq),
                getFastqRecords(fq1),
                getFastqRecords(fq2));
        Assert.assertEquals(result.get(1).size(), result.get(2).size());
        for (int i = 0; i < result.get(1).size(); i++) {
            Assert.assertEquals(result.get(1).get(i).getReadName(), result.get(2).get(i).getReadName());
        }
        return result;
    }
    @Test
    public void shouldFlipIfOnNegativeStrand() throws IOException {
        SAMRecord r1 = withName("r1", onNegative(withSequence("AACC", Read(0, 1, "4M"))))[0];
        List<List<FastqRecord>> result = go(ImmutableList.of("r1"), r1);
        Assert.assertEquals("GGTT", result.get(0).get(0).getReadString());
    }
    @Test
    public void shouldPairRecords() throws IOException {
        SAMRecord[] r1 = withName("r1", DP(0, 1, "10M", true, 1, 10, "10M", true));
        SAMRecord[] r2 = withName("r2", DP(0, 01, "10M", true, 0, 20, "10M", true));
        List<List<FastqRecord>> result = go(ImmutableList.of("r1", "r2"), r1[0], r2[0], r2[1], r1[1]);

        Assert.assertEquals("r2", result.get(1).get(0).getReadName());
        Assert.assertEquals("r1", result.get(1).get(1).getReadName());
    }
}