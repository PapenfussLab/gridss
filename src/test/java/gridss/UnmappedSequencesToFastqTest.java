package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class UnmappedSequencesToFastqTest extends IntermediateFilesTest {
    private List<FastqRecord> go(boolean includeSoftClips, int minLen, boolean uniqueNames, SAMRecord... records) {
        createInput(records);
        UnmappedSequencesToFastq cmd = new UnmappedSequencesToFastq();
        cmd.INPUT = input;
        cmd.OUTPUT = output;
        cmd.INCLUDE_SOFT_CLIPPED_BASES = includeSoftClips;
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
        List<FastqRecord> result = go(true, 1, false, r);
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
        List<FastqRecord> result = go(true, 1, false, r);
        Assert.assertEquals("GTTA", result.get(0).getReadString());
    }
    @Test
    public void hardClips_should_be_Ns() {
        SAMRecord r = Read(0, 1, "2H2S3M");
        List<FastqRecord> result = go(true, 1, false, r);
        Assert.assertEquals("NNAA", result.get(0).getReadString());
    }
    @Test
    public void shouldExportSoftClippedReads() {
        SAMRecord r1 = Read(0, 1, "5M5S");
        SAMRecord r2 = Read(0, 1, "6M4S");
        r1.setReadBases(B("CTTGGACGTA"));
        List<FastqRecord> result = go(true, 5, false, r1, r2);
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

        List<FastqRecord> result = go(true, 5, false, r1, r2, r3);
        Assert.assertEquals(2, result.size());
    }
    @Test
    public void shouldMakeUniqueName() {
        SAMRecord[] dp = DP(0, 1, "10M", true, 1, 1, "10M",true);
        for (SAMRecord r : dp) {
            r.setReadUnmappedFlag(true);
            r.setMappingQuality(0);
        }
        List<FastqRecord> result = go(true, 5, true, dp);
        Assert.assertNotEquals(result.get(0).getReadName(), result.get(1).getReadName());
    }
}