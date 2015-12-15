package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;


public class VcfBreakendToReadPairTest extends IntermediateFilesTest {
	private void create(ProcessingContext pc, Collection<? extends IdsvVariantContext> col) {
		VariantContextWriter vcw = pc.getVariantContextWriter(output, false);
		for (IdsvVariantContext c : col) {
			vcw.add(c);
		}
		vcw.close();
	}
	private File bam() { return new File(output.toString() + ".bam"); }
	private File bamFiltered() { return new File(output.toString() + ".filtered.bam"); }
	private VcfBreakendToReadPair getCmd() {
		VcfBreakendToReadPair cmd = new VcfBreakendToReadPair();
		cmd.INPUT = output;
		cmd.OUTPUT = bam();
		cmd.OUTPUT_FILTERED = bamFiltered();
		cmd.REFERENCE = getCommandlineContext().getReferenceFile();
		return cmd;
	}
	public List<List<SAMRecord>> go(ProcessingContext pc, Collection<? extends IdsvVariantContext> col) {
		create(pc, col);
		getCmd().doWork();
		return ImmutableList.of(getRecords(bam()), getRecords(bamFiltered()));
	}
	@Test
	public void should_write_sorted_read_per_breakend() throws IOException {
		ProcessingContext pc = getCommandlineContext();
		List<IdsvVariantContext> in = new ArrayList<IdsvVariantContext>();
		in.add(new IdsvVariantContextBuilder(pc, BP("ao", new BreakpointSummary(0, FWD, 1, 2, 3, BWD, 4, 5))).phredScore(6).make());
		in.add(new IdsvVariantContextBuilder(pc, BP("ah", new BreakpointSummary(3, BWD, 4, 5, 0, FWD, 1, 2))).phredScore(6).make());
		in.add(new IdsvVariantContextBuilder(pc, BP("bo", new BreakpointSummary(1, BWD, 10, 20, 2, FWD, 40, 50))).phredScore(100).make());
		in.add(new IdsvVariantContextBuilder(pc, BP("bh", new BreakpointSummary(2, FWD, 40, 50, 1, BWD, 10, 20))).phredScore(100).make());
		List<SAMRecord> list = go(pc, in).get(0);
		
		assertEquals(4, list.size());
		assertTrue(Ordering.from(new SAMRecordCoordinateComparator()).isOrdered(list));
		assertEquals(0, (int)list.get(0).getReferenceIndex());
		assertEquals(1, list.get(0).getAlignmentStart());
		assertEquals(2, list.get(0).getAlignmentEnd());
		assertFalse(list.get(0).getReadNegativeStrandFlag());
		assertTrue(list.get(2).getReadNegativeStrandFlag());
		assertTrue(list.get(0).getReadPairedFlag());
	}
	@Test
	public void should_match_breakends_using_event() throws IOException {
		ProcessingContext pc = getCommandlineContext();
		List<IdsvVariantContext> in = new ArrayList<IdsvVariantContext>();
		in.add(new IdsvVariantContextBuilder(pc, BP("test1o", new BreakpointSummary(0, FWD, 1, 2, 3, BWD, 4, 5)))
			.phredScore(6)
			.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, "test1")
			.make());
		in.add(new IdsvVariantContextBuilder(pc, BP("test1h", new BreakpointSummary(3, BWD, 4, 5, 0, FWD, 1, 2)))
			.phredScore(6)
			.attribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, "test1")
			.make());
		List<SAMRecord> list = go(pc, in).get(0);
		assertEquals("test1", list.get(0).getReadName());
		assertEquals("test1", list.get(1).getReadName());
		assertTrue(list.get(0).getFirstOfPairFlag());
		assertFalse(list.get(1).getFirstOfPairFlag());
		assertFalse(list.get(0).getSecondOfPairFlag());
		assertTrue(list.get(1).getSecondOfPairFlag());
	}
}
