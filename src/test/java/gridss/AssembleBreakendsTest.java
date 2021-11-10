package gridss;


import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import com.google.common.util.concurrent.MoreExecutors;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class AssembleBreakendsTest extends IntermediateFilesTest {
    @Test
    public void line_regression() throws IOException {
        File dir = new File("src/test/resources/anchor_misassembly/");
        AssembleBreakends cmd = new AssembleBreakends();
        cmd.OUTPUT = output;
        cmd.INPUT = ImmutableList.of(new File(dir, "anchor_misassembly.bam"));
        cmd.setReference(new File(dir, "ref.fa"));
        cmd.WORKING_DIR = dir;
        cmd.TMP_DIR = ImmutableList.of(dir);
        cmd.doWork(MoreExecutors.newDirectExecutorService());
        List<SAMRecord> asm = getRecords(output);
        Assert.assertEquals(2, asm.size());
        Assert.assertFalse(asm.get(1).getReadString().contains("AAAAAAAAAAAAAAAAAATGGTGGGACTCAGG"));
    }
}