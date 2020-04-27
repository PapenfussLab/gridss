package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.io.Files;
import gridss.cmdline.ReferenceCommandLineProgram;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

public class PrepareReferenceTest extends IntermediateFilesTest {
    @Test
    public void should_create_dictionary_bwaimage_gridsscache() throws IOException {
        File fa = testFolder.newFile("test.fa");
        File fai = testFolder.newFile("test.fa.fai");
        File dict = new File(testFolder.getRoot(), "test.fa.dict");
        File img = new File(testFolder.getRoot(), "test.fa.img");
        File cache = new File(testFolder.getRoot(), "test.fa.gridsscache");
        Files.copy(SMALL_FA_FILE, fa);
        Files.copy(new File(SMALL_FA_FILE.getAbsolutePath() + ".fai"), fai);
        PrepareReference cmd = new PrepareReference();
        cmd.instanceMain(new String[] {
                "R=" + fa.getAbsolutePath(),
        });
        assertTrue(dict.exists());
        assertTrue(img.exists());
        assertTrue(cache.exists());
        assertTrue(dict.length() > 0);
        assertTrue(img.length() > 0);
        assertTrue(cache.length() > 0);
    }
}