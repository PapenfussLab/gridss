package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class UnmappedSequencesToFastqTest extends IntermediateFilesTest {
    private List<FastqRecord> go(boolean includeSoftClips, int minLen, boolean uniqueNames, boolean internalBases, SAMRecord... records) {
        createInput(records);
        UnmappedSequencesToFastq cmd = new UnmappedSequencesToFastq();
        cmd.INPUT = ImmutableList.of(input);
        cmd.OUTPUT = output;
        cmd.INCLUDE_SOFT_CLIPPED_BASES = includeSoftClips;
        cmd.INCLUDE_UNMAPPED_INTERNAL_BASES = internalBases;
        cmd.MIN_SEQUENCE_LENGTH = minLen;
        cmd.UNIQUE_NAME = uniqueNames;
        cmd.doWork();
        Assert.assertTrue(output.exists());
        return getFastqRecords(output);
    }
    @Test
    public void shouldNotExportMappedSplitReads() {
        SAMRecord r = Read(0, 1, "50M50S");
        r.setAttribute("SA", "polyA,1,+,50S50M,0,0");
        List<FastqRecord> result = go(true, 1, false, false, r);
        Assert.assertEquals(0, result.size());
    }
    @Test
    public void shouldExportLargestUnmappedSequence() {
        //AACCGGTTACGT
        //SSMMM
        //         M
        SAMRecord r = Read(0, 1, "2S3M7S");
        r.setReadBases(B("AACCGGTTACGT"));
        r.setAttribute("SA", "polyA,1,+,9S1M2S,0,0");
        List<FastqRecord> result = go(true, 1, false, true, r);
        Assert.assertEquals("GTTA", result.get(0).getReadString());
    }
    @Test
    public void hardClips_should_be_Ns() {
        SAMRecord r = Read(0, 1, "2H2S3M");
        List<FastqRecord> result = go(true, 1, false, false, r);
        Assert.assertEquals("NNAA", result.get(0).getReadString());
    }
    @Test
    public void shouldExportSoftClippedReads() {
        SAMRecord r1 = Read(0, 1, "5M5S");
        SAMRecord r2 = Read(0, 1, "6M4S");
        r1.setReadBases(B("CTTGGACGTA"));
        List<FastqRecord> result = go(true, 5, false, false, r1, r2);
        Assert.assertEquals(1, result.size());
        Assert.assertEquals("ACGTA", result.get(0).getReadString());
    }
    @Test
    public void shouldExportUnmappedReads() {
        SAMRecord r1 = Read(0, 1, "10M");
        SAMRecord r2 = Read(0, 1, "10M");
        SAMRecord r3 = Read(1, 1, "10M");

        // TODO: do we care about strand?
        for (SAMRecord r : new SAMRecord[] { r1, r3}) {
            r.setReadUnmappedFlag(true);
            r.setMappingQuality(0);
        }

        List<FastqRecord> result = go(true, 5, false, false, r1, r2, r3);
        Assert.assertEquals(2, result.size());
    }
    @Test
    public void shouldMakeUniqueName() {
        SAMRecord[] dp = DP(0, 1, "10M", true, 1, 1, "10M",true);
        for (SAMRecord r : dp) {
            r.setReadUnmappedFlag(true);
            r.setMappingQuality(0);
        }
        List<FastqRecord> result = go(true, 5, true, false, dp);
        Assert.assertNotEquals(result.get(0).getReadName(), result.get(1).getReadName());
    }
    @Test
    public void shouldExtractQuals() {
        SAMRecord r = Read(0, 1, "1M4S");
        r.setBaseQualities(new byte[] { 1, 2, 3, 4, 5});
        List<FastqRecord> result = go(true, 1, false, false, r);
        Assert.assertArrayEquals(new byte[] { 2,3,4,5}, result.get(0).getBaseQualities());
    }
    @Test
    public void should_output_internal_sequence() throws IOException {
        UnmappedSequencesToFastq cmd = new UnmappedSequencesToFastq();
        cmd.INCLUDE_UNMAPPED_INTERNAL_BASES = true;
        file_test(cmd, new File("src/test/resources/unmappedSequencesToFastq/internal.sam")); }
    @Test
    public void should_match_c_implementation() throws IOException {
        UnmappedSequencesToFastq cmd = new UnmappedSequencesToFastq();
        file_test(cmd, new File("src/test/resources/unmappedSequencesToFastq/in.sam"));
    }
    private void file_test(UnmappedSequencesToFastq cmd, File in) throws IOException {
        cmd.INPUT = ImmutableList.of(in);
        cmd.OUTPUT = output;
        cmd.doWork();
        Assert.assertTrue(output.exists());
        assertLinesMatch(new File(in.toString() + ".out.fq"), output);
    }
}