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
        SAMRecord r = withSequence("NAAAAAAN", Read(0, 1, "1M7S"))[0];
        SAMRecord ra = withReadName("#1",withSequence("AAAAAAN", Read(0, 10, "7M")))[0];
        SplitReadHelper.convertToSplitRead(r, ImmutableList.of(ra), null, false);
        List<DirectedEvidence> evidence = ImmutableList.of(r, ra).stream()
            .map(x -> (DirectedEvidence)SingleReadEvidence.createEvidence(ses, 0, x).get(0))
            .sorted(DirectedEvidenceOrder.ByNatural)
            .collect(Collectors.toList());
        // Call variants
        List<VariantContextDirectedEvidence> calls = Lists.newArrayList(new VariantCallIterator(pc, evidence.iterator()));

        assertEquals(new BreakpointSummary(new BreakendSummary(0, FWD, 4, 1, 7), new BreakendSummary(0, BWD, 13, 10, 16)), calls.get(0).getBreakendSummary());
        assertEquals(calls.get(0).getBreakendSummary(), ((DirectedBreakpoint)calls.get(1)).getBreakendSummary().remoteBreakpoint());
    }
    @Test
    public void should_force_matching_nominal_positions() {
        ProcessingContext pc = getContext();
        MockSAMEvidenceSource ses = SES(pc);

        SAMRecord r1 = withSequence("NAAAAAAN", Read(0, 1, "1M7S"))[0];
        SAMRecord r1a = withReadName("#1",withSequence("AAAAAAN", Read(0, 10, "7M")))[0];
        SAMRecord r2 = withSequence("NAAAAAAN", Read(0, 1, "7S1M"))[0];
        SAMRecord r2a = withReadName("#1",withSequence("NAAAAAA", Read(0, 1, "7M")))[0];
        SplitReadHelper.convertToSplitRead(r1, ImmutableList.of(r1a), null, false);
        SplitReadHelper.convertToSplitRead(r2, ImmutableList.of(r2a), null, false);
        List<DirectedEvidence> evidence = ImmutableList.of(r1, r2).stream()
                .map(r -> (DirectedEvidence)SingleReadEvidence.createEvidence(ses, 0, r).get(0))
                .sorted(DirectedEvidenceOrder.ByNatural)
                .collect(Collectors.toList());
        // Set up calls with different nominal positions
        List<VariantContextDirectedEvidence> calls = ImmutableList.of(
                (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(pc).phredScore(1).attribute("EVENT", "event").breakpoint(new BreakpointSummary(0, FWD, 3, 3, 7, 0, BWD, 6, 6, 10), "").make(),
                (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(pc).phredScore(1).attribute("EVENT", "event").breakpoint(new BreakpointSummary(0, BWD, 10, 6, 10, 0, FWD, 7, 3, 7), "").make());

        ArrayList<VariantContextDirectedEvidence> annotatatedCalls = Lists.newArrayList(new SequentialEvidenceAnnotator(
                pc, calls.iterator(), evidence.iterator(), Collections.emptyIterator(), 1000, true, 1, MoreExecutors.newDirectExecutorService()));
        assertEquals(2, annotatatedCalls.size());
        assertEquals(annotatatedCalls.get(0).getBreakendSummary(), ((DirectedBreakpoint)annotatatedCalls.get(1)).getBreakendSummary().remoteBreakpoint());
    }
    @Test
    public void should_report_reference_base_at_nominal_position() {
        ProcessingContext pc = getContext();
        assertEquals("A", new IdsvVariantContextBuilder(pc).breakend(new BreakendSummary(1, FWD, 1), "").make().getReference().getBaseString());
        assertEquals("C", new IdsvVariantContextBuilder(pc).breakend(new BreakendSummary(1, FWD, 2, 1, 5), "").make().getReference().getBaseString());
        assertEquals("G", new IdsvVariantContextBuilder(pc).breakend(new BreakendSummary(1, FWD, 3, 1, 5), "").make().getReference().getBaseString());
        assertEquals("T", new IdsvVariantContextBuilder(pc).breakend(new BreakendSummary(1, FWD, 4, 1, 5), "").make().getReference().getBaseString());
    }
}