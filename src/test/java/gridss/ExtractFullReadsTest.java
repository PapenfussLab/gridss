package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.io.Files;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class ExtractFullReadsTest extends IntermediateFilesTest {
    public File simpleSetup() throws IOException {
        List<SAMRecord> in = new ArrayList<>();
        for (int i = 0 ; i < 1000; i++) {
            in.add(withName(String.format("a%d", i), Read(0, 1, "100M"))[0]); // overlaps
            in.add(withName(String.format("b%d", i), Read(0, 100, "100M"))[0]); // off by 1
            in.add(withName(String.format("c%d", i), Read(0, 150, "100M"))[0]);
            in.add(withName(String.format("d%d", i), Read(1, 100, "100M"))[0]);
            in.add(withName(String.format("e%d", i), Read(1, 200, "100M"))[0]); // in by 1
            in.add(withName(String.format("f%d", i), Read(1, 201, "100M"))[0]); // overlaps
        }
        createInput(in);

        File bed = new File(testFolder.getRoot(), "region.bed");
        try(BufferedWriter writer = Files.newWriter(bed, StandardCharsets.UTF_8)) {
            writer.write("polyA 0   100\n");
            writer.write("polyACGT 298   1000\n");
        }
        return bed;
    }
    @Test
    public void test_linear_scan() throws IOException {
        File bed = simpleSetup();
        try(BufferedWriter writer = Files.newWriter(bed, StandardCharsets.UTF_8)) {
            writer.write("polyA 0   99\n");
            writer.write("polyACGT 201   1000\n");
        }
        String[] args = new String[] {
                "INPUT=" + input.toString(),
                "REFERENCE_SEQUENCE=" + reference.toString(),
                "OUTPUT=" + output.toString(),
                "REGION_BED=" + bed.toString()
        };
        assertEquals(0, new ExtractFullReads().instanceMain(args));
        bed.delete();
        List<SAMRecord> outputRecords = getRecords(output);
        assertEquals(3000, outputRecords.size());
    }
}