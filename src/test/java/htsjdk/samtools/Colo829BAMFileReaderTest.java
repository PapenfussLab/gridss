package htsjdk.samtools;

import htsjdk.samtools.util.CloseableIterator;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

public class Colo829BAMFileReaderTest {
    private static File TEST_FILE = new File("D:/dev/colo829/chr1/colo829n_chr1.bam");
    private static File TEST_FILE_INDEX = new File("D:/dev/colo829/chr1/colo829n_chr1.bai");
    @Test
    public void index_read_should_match_read_counts() throws IOException {
        QueryInterval[] qi = new QueryInterval[1024];
        for (int i = 0; i < 1024; i++) {
            int start = i * (249200621 / 1024);
            qi[i] = new QueryInterval(0, start, start + 1000);
        }
        BAMFileReader sync = new BAMFileReader(TEST_FILE, TEST_FILE_INDEX, false, false, ValidationStringency.SILENT, new DefaultSAMRecordFactory());
        CloseableIterator<SAMRecord> it = sync.query(qi, false);
        int syncCount = 0;
        while (it.hasNext()) {
            it.next();
            syncCount++;
        }
        BAMFileReader async = new BAMFileReader(TEST_FILE, TEST_FILE_INDEX, false, true, ValidationStringency.SILENT, new DefaultSAMRecordFactory());
        CloseableIterator<SAMRecord> ait = async.query(qi, false);
        int asyncCount = 0;
        while (it.hasNext()) {
            it.next();
            syncCount++;
        }
        sync.close();
        async.close();
        Assert.assertEquals(syncCount, asyncCount);
    }
}