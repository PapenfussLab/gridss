package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;

public class SequentialEvidenceAnnotatorTest extends TestHelper {
    private static List<SAMRecord> getSplitReadTestData(ProcessingContext pc) {
        // 12345678901234567890
        //    AAAA
        //      AAAA
        // nnn      nnn
        String seq = "nnnAAAAnnn";
        SAMRecord[] left = new SAMRecord[] {
                withName("left", withSequence(seq, Read(0, 1, "3M7S")))[0],
                withName("left#0", withSequence(seq, Read(0, 1, "3S7M")))[0]
        };
        SAMRecord[] right = new SAMRecord[] {
                withMapq(40, withName("right", withSequence(seq, Read(0, 1, "7M3S"))))[0],
                withMapq(40, withName("right#0", withSequence(seq, Read(0, 1, "7S3M"))))[0]
        };
        SplitReadHelper.convertToSplitRead(left[0], ImmutableList.of(left[1]), pc.getReference(), false);
        SplitReadHelper.convertToSplitRead(right[0], ImmutableList.of(right[1]), pc.getReference(), false);
        List<SAMRecord> reads = ImmutableList.of(left[0], left[1], right[0], right[1]);
        return reads;
    }
    @Test
    public void should_centre_align() {
        ProcessingContext pc = getContext();
        MockSAMEvidenceSource ses = SES(pc);
        List<SAMRecord> reads = getSplitReadTestData(pc);
        List<DirectedEvidence> evidence = reads.stream()
            .map(r -> (DirectedEvidence)SingleReadEvidence.createEvidence(ses, 0, r).get(0))
            .sorted(DirectedEvidenceOrder.ByNatural)
            .collect(Collectors.toList());
        // Call variants
        List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new VariantCallIterator(pc, evidence.iterator()));
        assertEquals(new BreakpointSummary(new BreakendSummary(0, FWD, 5, 5, 7), new BreakendSummary(0, BWD, 8, 6, 10)), calls.get(0).getBreakendReadCount());
        assertEquals(calls.get(0).getBreakendSummary(), ((DirectedBreakpoint)calls.get(1)).getBreakendSummary().remoteBreakpoint());
    }
    @Test
    public void should_force_matching_nominal_positions() {
        ProcessingContext pc = getContext();
        MockSAMEvidenceSource ses = SES(pc);
        List<SAMRecord> reads = getSplitReadTestData(pc);
        List<DirectedEvidence> evidence = reads.stream()
                .map(r -> (DirectedEvidence)SingleReadEvidence.createEvidence(ses, 0, r).get(0))
                .sorted(DirectedEvidenceOrder.ByNatural)
                .collect(Collectors.toList());
        // Set up calls with different nominal positions
        List<VariantContextDirectedEvidence> calls = ImmutableList.of(
                (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(pc).breakpoint(new BreakpointSummary(0, FWD, 3, 3, 7, 0, BWD, 6, 6, 10), "").make(),
                (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(pc).breakpoint(new BreakpointSummary(0, BWD, 6, 10, 10, 0, FWD, 3, 7, 7), "").make());

        ArrayList<VariantContextDirectedEvidence> annotatatedCalls = Lists.newArrayList(new SequentialEvidenceAnnotator(
                pc, calls.iterator(), Collections.emptyIterator(), Collections.emptyIterator(), 1000, true, 0, MoreExecutors.newDirectExecutorService()));
        assertEquals(2, annotatatedCalls.size());
        assertEquals(calls.get(0).getBreakendSummary(), ((DirectedBreakpoint)calls.get(1)).getBreakendSummary().remoteBreakpoint());
    }
    @Test
    public void should_report_reference_base_at_nominal_position() {
        ProcessingContext pc = getContext();
        MockSAMEvidenceSource ses = SES(pc);
        // 12345678901234567890
        // ACGTACGTACGTACGT
        // nnnn    dddd
        //     MMMM    NNNN
        String seq = "nnnnACGTACGTnnnn";
        List<SAMRecord> reads = Lists.newArrayList(withSequence(seq, Read(0, 1, "8M4D4M")));
        List<DirectedEvidence> evidence = reads.stream()
                .map(r -> (DirectedEvidence) SingleReadEvidence.createEvidence(ses, 0, r).get(0))
                .sorted(DirectedEvidenceOrder.ByNatural)
                .collect(Collectors.toList());
        List<VariantContextDirectedEvidence> calls = ImmutableList.of(
                (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(pc).breakpoint(new BreakpointSummary(0, FWD, 5,  0, BWD, 10), "").make());
        assertEquals("A", calls.get(0).getReference().getBaseString());

        ArrayList<VariantContextDirectedEvidence> annotatatedCalls = Lists.newArrayList(new SequentialEvidenceAnnotator(
                pc, calls.iterator(), evidence.iterator(), Collections.emptyIterator(), 1000, true, 0, MoreExecutors.newDirectExecutorService()));
        // Sanity precondition that we're centre-aligned
        assertEquals(annotatatedCalls.get(0).getBreakendSummary(), new BreakpointSummary(0, FWD, 6, 4, 8, 0, BWD, 11, 9, 13));
        // check we are reporting the centre-aligned reference base
        assertEquals("C", annotatatedCalls.get(0).getReference().getBaseString());
    }
}