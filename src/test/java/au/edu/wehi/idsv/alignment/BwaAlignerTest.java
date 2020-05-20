package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BwaAlignerTest extends TestHelper {
    public static List<SAMRecord> hitsFor(String readname, List<SAMRecord> hits) {
        return hits.stream().filter(x -> x.getReadName().equals(readname)).collect(Collectors.toList());
    }
    @Test
    public void should_create_index_from_reference() throws IOException {
        TemporaryFolder testFolder = new TemporaryFolder();
        try {
            testFolder.create();
            File ref = new File(testFolder.getRoot(), "small.fa");
            File image = new File(testFolder.getRoot(), "small.fa.img");
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa"), testFolder.getRoot());
            BwaAligner.createBwaIndexFor(ref);
            assertTrue(image.exists());
        } finally {
            testFolder.delete();
        }
    }
    @Test
    public void should_lazily_create_from_bwa_index_files() throws IOException {
        TemporaryFolder testFolder = new TemporaryFolder();
        try {
            testFolder.create();
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa"), testFolder.getRoot());
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa.amb"), testFolder.getRoot());
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa.ann"), testFolder.getRoot());
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa.bwt"), testFolder.getRoot());
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa.pac"), testFolder.getRoot());
            FileUtils.copyFileToDirectory(new File("src/test/resources/small.fa.sa"), testFolder.getRoot());
            File ref = new File(testFolder.getRoot(), "small.fa");
            File image = new File(testFolder.getRoot(), "small.fa.img");
            BwaAligner.createBwaIndexFor(ref);
            assertTrue(image.exists());
        } finally {
            testFolder.delete();
        }
    }
    @Test(expected = RuntimeException.class)
    public void should_except_if_missing_index() {
        try (BwaAligner ba = new BwaAligner(new File("src/test/resources/missing_index"), SMALL_FA.getSequenceDictionary(), 2)) {
            List<SAMRecord> result = ba.align(ImmutableList.of(new FastqRecord("read", "A", "", "A")));
        }
    }
    @Test
    public void should_run_bwa_jni() throws IOException {
        try (BwaAligner ba = new BwaAligner(new File("src/test/resources/small.fa"), SMALL_FA.getSequenceDictionary(), 2)) {
            Collection<FastqRecord> input = ImmutableList.of(
                    new FastqRecord("noHit", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
                    new FastqRecord("polyA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                    new FastqRecord("polyACGT", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
                    new FastqRecord("random", S(RANDOM).substring(0, 100), "", S(getPolyA(100))));
            List<SAMRecord> result = ba.align(input);
            assertEquals(1, hitsFor("noHit", result).size());
            assertTrue(hitsFor("noHit", result).get(0).getReadUnmappedFlag());

            assertEquals(1, hitsFor("random", result).size());
            assertEquals("100M", hitsFor("random", result).get(0).getCigarString());
            assertEquals("random", hitsFor("random", result).get(0).getReferenceName());
            assertEquals(1, hitsFor("random", result).get(0).getAlignmentStart());
        }
    }
    @Test
    public void should_return_soft_clips() throws IOException {
        try (BwaAligner ba = new BwaAligner(new File("src/test/resources/small.fa"), SMALL_FA.getSequenceDictionary(), 2)) {
            Collection<FastqRecord> input = ImmutableList.of(
                    new FastqRecord("random", "CCCCCCCCCCCCCCCCCCCC" + S(RANDOM).substring(0, 100), "", S(getPolyA(120))));
            List<SAMRecord> result = ba.align(input);
            assertEquals(1, result.size());
            assertEquals("20S100M", result.get(0).getCigarString());
            assertEquals(1, result.get(0).getAlignmentStart());
        }
    }
    @Test
    public void should_have_min_alignment_length() throws IOException {
        String s = S(getPolyA(120));
        try (BwaAligner ba = new BwaAligner(new File("src/test/resources/small.fa"), SMALL_FA.getSequenceDictionary(), 2)) {
            List<FastqRecord> input = new ArrayList<>();
            for (int i = 0; i < 100; i++) {
                input.add(new FastqRecord("read" + i, s.substring(0, i), "", s.substring(0, i)));
            }
            List<SAMRecord> result = ba.align(input);
            for (int i = 0; i < 100; i++) {
                Assert.assertEquals(i, result.get(i).getReadLength());
            }
        }
    }
    @Test
    public void should_handle_Ns() throws IOException {
        String s = S(getPolyA(120));
        try (BwaAligner ba = new BwaAligner(new File("src/test/resources/small.fa"), SMALL_FA.getSequenceDictionary(), 2)) {
            List<FastqRecord> input = new ArrayList<>();
            for (int i = 0; i < 100; i++) {
                input.add(new FastqRecord("read" + i, getPoly(i, 'N'), "", getPolyA(i)));
            }
            List<SAMRecord> result = ba.align(input);
            for (int i = 0; i < 100; i++) {
                Assert.assertEquals(i, result.get(i).getReadLength());
            }
        }
    }
}