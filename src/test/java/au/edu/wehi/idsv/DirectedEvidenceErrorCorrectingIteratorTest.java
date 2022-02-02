package au.edu.wehi.idsv;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class DirectedEvidenceErrorCorrectingIteratorTest extends TestHelper {
    @Test
    public void should_not_EC_line_anchor() throws IOException {
        File dir = new File("src/test/resources/anchor_misassembly/");
        File ref = new File(dir, "ref.fa");
        ReferenceLookup rl = new TwoBitBufferedReferenceSequenceFile(new IndexedFastaSequenceFile(ref));
        FileSystemContext fsc = new FileSystemContext(dir, 500000);
        ProcessingContext pc = new ProcessingContext(fsc, ref, rl, new ArrayList<Header>(), getConfig());
        pc.registerCategory("test");
        SAMEvidenceSource ses = new SAMEvidenceSource(pc, new File(dir, "anchor_misassembly.bam"), null, 0);
        DirectedEvidenceErrorCorrectingIterator deeci = new DirectedEvidenceErrorCorrectingIterator(
                pc.getLinear(),
                ses.iterator(SAMEvidenceSource.EvidenceSortOrder.SAMRecordStartPosition),
                ses.getMinConcordantFragmentSize(),
                ses.getMaxConcordantFragmentSize(),
                ses.getMaxReadLength(),
                ses.getMaxReadMappedLength(),
                SAMEvidenceSource.EvidenceSortOrder.SAMRecordStartPosition,
                10,
                21,
                100,
                false);
        List<DirectedEvidence> de = Lists.newArrayList(deeci);
        List<SAMRecord> inSam = ses.iterator(SAMEvidenceSource.EvidenceSortOrder.SAMRecordStartPosition).stream().map(e -> {
                try {
                    return (SAMRecord)e.getUnderlyingSAMRecord().clone();
                } catch (CloneNotSupportedException ex) {
                    throw new RuntimeException(ex);
                }
            }).collect(Collectors.toList());
        List<SAMRecord> outSam = de.stream()
                .map(e -> e.getUnderlyingSAMRecord())
                .collect(Collectors.toList());
        // check that some error correction has been done
        Assert.assertFalse(Streams.zip(inSam.stream(), outSam.stream(), (a, b) -> Pair.of(a, b))
                .allMatch(pair -> pair.getLeft().equals(pair.getRight())));
        Assert.assertTrue(outSam.stream()
                .filter(s -> s.getAlignmentStart() == 1081)
                .allMatch(s -> new String(s.getReadBases()).contains("AAGAACAAGTGGTGGGACT")));
    }
}