package au.edu.wehi.idsv.pipeline;

import static org.junit.Assert.assertEquals;

import java.util.EnumSet;
import java.util.List;

import org.junit.Test;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;


public class SortRealignedSoftClipsTest extends IntermediateFilesTest {
	ProcessingContext processContext;
	SAMEvidenceSource source;
	public void go(boolean perChr, SAMRecord... realign) {
		processContext = getCommandlineContext(perChr);
		processContext.getRealignmentParameters().minLength = 15;
		source = new SAMEvidenceSource(processContext, input, 0);
		source.completeSteps(ProcessStep.ALL_STEPS);
		if (perChr) {
			for (final SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				List<SAMRecord> forchr = Lists.newArrayList(Iterables.filter(Lists.newArrayList(realign), new Predicate<SAMRecord>() {
					@Override
					public boolean apply(SAMRecord arg0) {
						return Integer.parseInt(arg0.getReadName().split("#")[0]) == seq.getSequenceIndex();
					}
				}));
				createBAM(processContext.getFileSystemContext().getRealignmentBamForChr(input, seq.getSequenceName(), 0), SortOrder.unsorted, forchr.toArray(new SAMRecord[forchr.size()]));
			}
		} else {
			createBAM(processContext.getFileSystemContext().getRealignmentBam(input, 0), SortOrder.unsorted, realign);
		}
		SortRealignedSoftClips srs = new SortRealignedSoftClips(processContext, source);
		srs.process(EnumSet.allOf(ProcessStep.class));
		srs.close();
	}
	@Test
	public void should_sort_by_realignment_position() {
		createInput(
				withReadName("r1", Read(0, 1, "15M15S")),
				withReadName("r2", Read(1, 2, "15M15S")),
				withReadName("r3", Read(2, 3, "15M15S")));
		go(false,
				withReadName("0#1#0#fr1", Read(2, 10, "15M"))[0],
				withReadName("1#2#0#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#0#fr3", Read(0, 10, "15M"))[0]
		);
		assertEquals(3, getRSC(source).size());
		assertEquals("r3", getRSC(source).get(0).getReadName());
		assertEquals("r2", getRSC(source).get(1).getReadName());
		assertEquals("r1", getRSC(source).get(2).getReadName());
		
		assertEquals(3, getRRR(source).size());
		assertEquals("2#3#0#fr3", getRRR(source).get(0).getReadName());
		assertEquals("1#2#0#fr2", getRRR(source).get(1).getReadName());
		assertEquals("0#1#0#fr1", getRRR(source).get(2).getReadName());
	}
	@Test
	public void should_write_sc_for_both_forward_and_reverse_realignment() {
		SAMRecord[] realigned = new SAMRecord[] {
			withReadName("0#1#0#br1", Unmapped(15))[0],
			withReadName("0#1#0#fr1", Read(2, 10, "15M"))[0],
			withReadName("1#2#0#fr2", Read(1, 10, "15M"))[0],
			withReadName("1#2#0#br2", Read(1, 15, "15M"))[0],
			withReadName("2#3#0#br3", Unmapped(15))[0],
			withReadName("2#3#0#fr3", Read(0, 10, "15M"))[0],
		};
		SAMRecord[] softClip = new SAMRecord[] {
			withReadName("r1", Read(0, 1, "15S15M15S"))[0],
			withReadName("r2", Read(1, 2, "15S15M15S"))[0],
			withReadName("r3", Read(2, 3, "15S15M15S"))[0]
		};
		createInput(softClip);
		ProcessingContext processContext = getCommandlineContext(false);
		processContext.getRealignmentParameters().minLength = 1;
		SAMEvidenceSource source = new SAMEvidenceSource(processContext, input, 0);
		source.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(processContext.getFileSystemContext().getRealignmentBam(input, 0), SortOrder.unsorted, realigned);
		SortRealignedSoftClips srs = new SortRealignedSoftClips(processContext, source);
		srs.process(EnumSet.allOf(ProcessStep.class));
		srs.close();
		assertEquals(4, getRSC(source).size());
	}
	@Test
	public void should_sort_by_realignment_position_per_chr() {
		createInput(
				withReadName("r1", Read(0, 1, "15M15S")),
				withReadName("r2", Read(1, 2, "15M15S")),
				withReadName("r3", Read(2, 3, "15M15S")));
		go(true,
				withReadName("0#1#0#fr1", Read(2, 10, "15M"))[0],
				withReadName("1#2#0#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#0#fr3", Read(0, 10, "15M"))[0]
		);
		assertEquals(3, new PerChr(processContext).getRSC(source).size());
		assertEquals("r3", new PerChr(processContext).getRSC(source).get(0).getReadName());
		assertEquals("r2", new PerChr(processContext).getRSC(source).get(1).getReadName());
		assertEquals("r1", new PerChr(processContext).getRSC(source).get(2).getReadName());
		
		assertEquals(3, new PerChr(processContext).getRRR(source).size());
		assertEquals("2#3#0#fr3", new PerChr(processContext).getRRR(source).get(0).getReadName());
		assertEquals("1#2#0#fr2", new PerChr(processContext).getRRR(source).get(1).getReadName());
		assertEquals("0#1#0#fr1", new PerChr(processContext).getRRR(source).get(2).getReadName());
	}
	@Test
	public void should_not_include_spanned_indels() {
		createInput(withReadName("r1", Read(0, 1, "5M5D5M")));
		go(true);
		assertEquals(0, new PerChr(processContext).getRSC(source).size());
	}
}

