package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.IOException;
import java.util.function.Supplier;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BwaStreamingAlignerTest extends TestHelper {
    @Test
    @Category(JniAlignerTests.class)
    public void should_align_reads() throws IOException {
        BwaStreamingAligner bwamem = new BwaStreamingAligner(SMALL_FA_FILE, SMALL_FA.getSequenceDictionary(), 2, 1000);
        for (FastqRecord fq : ImmutableList.of(
                new FastqRecord("noHit", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
                new FastqRecord("polyA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                new FastqRecord("polyACGT", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
                new FastqRecord("random", S(RANDOM).substring(0, 100), "", S(getPolyA(100))))) {
            bwamem.asyncAlign(fq);
        }
        bwamem.flush();
        int i = 0;
        while (bwamem.processedAlignmentRecords() > 0) {
            bwamem.getAlignment();
            i++;
        }
        assertTrue(i >= 4);
    }
    @Test(expected = IllegalStateException.class)
    @Category(JniAlignerTests.class)
    public void should_batch_records_for_processing() throws IOException {
        BwaStreamingAligner bwamem = new BwaStreamingAligner(SMALL_FA_FILE, SMALL_FA.getSequenceDictionary(), 2, 1000);
        for (FastqRecord fq : ImmutableList.of(
                new FastqRecord("noHit", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
                new FastqRecord("polyA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                new FastqRecord("polyACGT", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
                new FastqRecord("random", S(RANDOM).substring(0, 100), "", S(getPolyA(100))))) {
            bwamem.asyncAlign(fq);
        }
        assertEquals(0, bwamem.processedAlignmentRecords());
        bwamem.getAlignment();
    }
    private static void waitUntilTrue(Supplier<Boolean> condition) throws InterruptedException {
        for (int i = 0; i < 1024; i++) {
            if (condition.get()) {
                break;
            } else {
                Thread.sleep(1);
            }
        }
    }
    @Test
    @Category(JniAlignerTests.class)
    public void counts_should_match_aligner_state_when_batch_processed() throws IOException, InterruptedException {
        final BwaStreamingAligner bwamem = new BwaStreamingAligner(SMALL_FA_FILE, SMALL_FA.getSequenceDictionary(), 2, 1);
        bwamem.asyncAlign(new FastqRecord("noHit", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
        waitUntilTrue(() -> bwamem.outstandingAlignmentRecord() == 0);
        assertEquals(0, bwamem.outstandingAlignmentRecord());
        assertEquals(1, bwamem.processedAlignmentRecords());
        bwamem.close();
    }
    @Test
    @Category(JniAlignerTests.class)
    public void flush_should_force_processing_and_block_till_completed() throws IOException, InterruptedException {
        final BwaStreamingAligner bwamem = new BwaStreamingAligner(SMALL_FA_FILE, SMALL_FA.getSequenceDictionary(), 2, 1000);
        bwamem.asyncAlign(new FastqRecord("noHit", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
        assertEquals(1, bwamem.outstandingAlignmentRecord());
        assertEquals(0, bwamem.processedAlignmentRecords());
        bwamem.flush();
        assertEquals(0, bwamem.outstandingAlignmentRecord());
        assertEquals(1, bwamem.processedAlignmentRecords());
        bwamem.getAlignment();
        assertEquals(0, bwamem.outstandingAlignmentRecord());
        assertEquals(0, bwamem.processedAlignmentRecords());
    }
}