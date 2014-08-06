package au.edu.wehi.idsv;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;


public class SortRealignedSoftClipsTest extends IntermediateFilesTest {
	ProcessingContext processContext;
	SAMEvidenceSource source;
	public void go(SAMRecord... realign) {
		processContext = getCommandlineContext(false);
		source = new SAMEvidenceSource(processContext, input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(processContext.getFileSystemContext().getRealignmentBam(input), SortOrder.unsorted, realign);
		SortRealignedSoftClips srs = new SortRealignedSoftClips(getCommandlineContext(), source);
		srs.process(false);
		
	}
	@Test
	public void should_sort_by_realignment_position() {
		createInput(
				withReadName("r1", Read(0, 1, "15M15S")),
				withReadName("r2", Read(1, 2, "15M15S")),
				withReadName("r3", Read(2, 3, "15M15S")));
		go(
				withReadName("0#1#fr1", Read(2, 10, "15M"))[0],
				withReadName("1#2#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#fr3", Read(0, 10, "15M"))[0]
		);
		assertEquals(3, getRSC(source).size());
		assertEquals("r3", getRSC(source).get(0).getReadName());
		assertEquals("r2", getRSC(source).get(1).getReadName());
		assertEquals("r1", getRSC(source).get(2).getReadName());
		
		assertEquals(3, getRRR(source).size());
		assertEquals("2#3#fr3", getRRR(source).get(0).getReadName());
		assertEquals("1#2#fr2", getRRR(source).get(1).getReadName());
		assertEquals("0#1#fr1", getRRR(source).get(2).getReadName());
	}
	@Test
	public void should_write_sc_for_both_forward_and_reverse_realignment() {
		SAMRecord[] realigned = new SAMRecord[] {
			withReadName("0#1#fr1", Read(2, 10, "15M"))[0],
			withReadName("1#2#fr2", Read(1, 10, "15M"))[0],
			withReadName("1#2#br2", Read(1, 15, "15M"))[0],
			withReadName("2#3#fr3", Read(0, 10, "15M"))[0],
		};
		SAMRecord[] softClip = new SAMRecord[] {
			withReadName("r1", Read(0, 1, "15S15M15S"))[0],
			withReadName("r2", Read(1, 2, "15S15M15S"))[0],
			withReadName("r3", Read(2, 3, "15S15M15S"))[0]
		};
		createInput(softClip);
		ProcessingContext processContext = getCommandlineContext(false);
		SAMEvidenceSource source = new SAMEvidenceSource(processContext, input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(processContext.getFileSystemContext().getRealignmentBam(input), SortOrder.unsorted, realigned);
		SortRealignedSoftClips srs = new SortRealignedSoftClips(processContext, source);
		srs.process(false);
		
		assertEquals(4, getRSC(source).size());
	}
}

