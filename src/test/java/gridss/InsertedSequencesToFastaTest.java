package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class InsertedSequencesToFastaTest extends IntermediateFilesTest {
    public static final File colo829somatic = new File("src/test/resources/gridss294.colo829.somatic.vcf");
    @Test
    public void shouldUseIDasName() {
        InsertedSequencesToFasta program = new InsertedSequencesToFasta();
        File output = new File(testFolder.getRoot(), "out.fa");
        program.INPUT = colo829somatic;
        program.OUTPUT = output;
        program.MIN_SEQUENCE_LENGTH = 10;
        program.doWork();

        Assert.assertTrue(output.exists());
        ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(output);
        Map<String, String> seq = new HashMap<>();
        for (ReferenceSequence s = rsf.nextSequence(); s != null; s = rsf.nextSequence()) {
            seq.put(s.getName(), s.getBaseString());
        }
        Assert.assertEquals("TTGTTGTTGTTGGTTTTTC", seq.get("gridss0fb_12021o"));
        Assert.assertEquals("TTGTTGTTGTTGGTTTTTC", seq.get("gridss0fb_12021h"));
        Assert.assertEquals("GAATCATCGAATGGTCTCGATTGGAATCATTATCAAATGGAATCGAATGGAATCACCGAATAGAATCGAATGGAACAATCATCGAATGGACTCAAATGGAATTATCCTCAAATGGAATGGAATGGAATTATTGAATGCAATCGAAAGGAATTATCGAATGCAATCGAATAGAATCATCGAATGGACTCGAATGGAATCATCGAATGGAATGGAATGGAACAGTCAATGAACTCGAATGGAATCATCATTGAATGGAATCGAATGGAATCATCGAGTGGAATCGAATGGAATTATGATCAAATGGAATCGAATGTAATCATCATCAAATGGAATCAAAAATAACCATCATCAATTGGTATTGAATGGAATTGTCATCAAATGGAATTCAAAGGAATCATCATCAAATGGAACCGAATGGAATCCTCATTGAATGGAAATGAAAGGGGTCATCATCTAATGGAATCGCATGGAATCATCATCAAATGGAATCGAATGTAATCTTCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGTAAA", seq.get("gridss51b_73915b"));
        Assert.assertFalse(seq.containsKey("gridss174ff_186h"));
    }
}