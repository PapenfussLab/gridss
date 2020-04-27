package au.edu.wehi.idsv;

import au.edu.wehi.idsv.alignment.*;
import au.edu.wehi.idsv.picard.BufferedReferenceSequenceFile;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class StreamingSplitReadRealignerTest extends SplitReadRealignerTest {
    private static final SmithWatermanStreamingAligner aligner = new SmithWatermanStreamingAligner(AlignerFactory.create(), SMALL_FA, 2);

    @Override
    protected SplitReadRealigner createAligner() {
        return new StreamingSplitReadRealigner(getContext(), aligner, 10);
    }

    public class StubStreamingAligner implements StreamingAligner {
        private int out = 0;
        private int in = 0;
        private SAMRecord[] alignments;
        public StubStreamingAligner(SAMRecord... alignments) {
            this.alignments = alignments;
        }
        @Override
        public void asyncAlign(FastqRecord fq) throws IOException {
            this.alignments[in++].setReadName(fq.getReadName());
        }
        @Override
        public void flush() throws IOException { }

        public boolean hasCompletedAlignmentRecord() {
            return in > out & out < alignments.length;
        }

        @Override
        public SAMRecord getAlignment() {
            return alignments[out++];
        }

        @Override
        public void close() throws IOException { }
        @Override
        public int processedAlignmentRecords() {
            // alignment is immediate
            if (!hasCompletedAlignmentRecord()) return 0;
            int records = in - out;
            return records;
        }
        @Override
        public int outstandingAlignmentRecord() {
            return 0;
        }
    }

    @Test
    public void streaming_should_realign_multiple_times() throws IOException {
        SAMRecord r0 = Read(0, 100, "10M30S");
        r0.setReadName("r");
        SAMRecord r1 = Read(0, 200, "10M20S");
        r1.setAttribute("SA", "this record is a split read realignment");
        SAMRecord r2 = Read(0, 200, "10M10S");
        SAMRecord r3 = Read(0, 200, "10M");
        createBAM(input, SAMFileHeader.SortOrder.coordinate, r0);
        StreamingSplitReadRealigner srr = new StreamingSplitReadRealigner(getContext(), new StubStreamingAligner(r1, r2, r3), 10);
        srr.createSupplementaryAlignments(input, output, output);
        List<SAMRecord> list = getRecords(output);
        assertEquals(4, list.size());
    }

    @Test
    @Category(ExternalAlignerTests.class)
    @Ignore("Working 2018-04-08. Currently need to manual check # restarts of external aligner to actually test this functionality. Needs a delayed return stub to test properly.")
    public void streaming_should_limit_outstanding_records_to_buffer_size() throws IOException, CloneNotSupportedException {
        ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
        BufferedReferenceSequenceFile lookup = new BufferedReferenceSequenceFile(ReferenceSequenceFileFactory.getReferenceSequenceFile(ExternalAlignerTests.REFERENCE));
        ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), ExternalAlignerTests.REFERENCE, lookup, Lists.newArrayList(), getConfig(testFolder.getRoot()));

        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(lookup.getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);


        SAMRecord r = new SAMRecord(header);
        r.setReferenceIndex(0);
        r.setAlignmentStart(1399998);
        r.setCigarString("301S103M");
        r.setReadBases(B("GGATATATAGGGATAGAAGCTTGAATAGTCTGGACATATATTTGTATTGAAATACAAATGTAAGATTTCAGTTAATCAATTTAAACATTTTTATTTTCAAGGGCTTCCAGCGTCCACTTCCTACGGCAAGCAGGAGGAGACAAGCGCCACCCTGCGCTCGCGGAGCCGACCCCGGCTCTCCCCTCCCGTGGCCGCAGGGGTCTGACAGAAAGGGGTCACTAATCTACTTGGCCTTTTGAGGACTGATCCTTAAGAATAATTTTTTTTTTTTTATGATCTTGAAGGCTGAGAAGTATTAGAGTAGGTTTTTTTCTCCTTCATAAGGCCAGATTCTTCTTTCTGTCACAGATTTCAAGTCCCCGCCTCAGCAGCCTTTCACTGTCAGTTCTTTCTCACGTGACCCT"));
        r.setBaseQualities(B("?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA"));
        r.setReadName("four_way_split_read");

        SAMRecord[] arr = new SAMRecord[64];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = (SAMRecord)r.clone();
            arr[i].setReadName(Integer.toString(i));
        }
        createBAM(input, header, arr);

        srr = new StreamingSplitReadRealigner(pc, aligner, 16);
        srr.createSupplementaryAlignments(input, output, output);
        List<SAMRecord> list = getRecords(output);
        assertEquals(4, list.size());
    }
}