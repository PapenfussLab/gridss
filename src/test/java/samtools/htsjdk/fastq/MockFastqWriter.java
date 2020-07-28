package samtools.htsjdk.fastq;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

import java.io.Closeable;
import java.util.ArrayList;
import java.util.List;

public class MockFastqWriter implements FastqWriter, Closeable {
    public List<FastqRecord> records = new ArrayList<>();
    @Override
    public void write(FastqRecord fastqRecord) {
        records.add(fastqRecord);
    }

    @Override
    public void close() {
    }
}
