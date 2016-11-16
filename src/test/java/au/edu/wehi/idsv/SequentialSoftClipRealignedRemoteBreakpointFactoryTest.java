package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.EnumSet;

import org.junit.Test;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.pipeline.SortRealignedSoftClips;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;


public class SequentialSoftClipRealignedRemoteBreakpointFactoryTest extends IntermediateFilesTest {
	@Test
	public void should_get_soft_clip_for_realignment() {
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
		
		SequentialSoftClipRealignedRemoteBreakpointFactory ffactory = new SequentialSoftClipRealignedRemoteBreakpointFactory(Iterators.peekingIterator(getRSC(source).iterator()), FWD);
		SequentialSoftClipRealignedRemoteBreakpointFactory bfactory = new SequentialSoftClipRealignedRemoteBreakpointFactory(Iterators.peekingIterator(getRSC(source).iterator()), BWD);
		
		SAMRecord[] orderedRealigned = new SAMRecord[] {
				withReadName("2#3#0#fr3", Read(0, 10, "15M"))[0],
				withReadName("1#2#0#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#0#br2", Read(1, 15, "15M"))[0],
				withReadName("0#1#0#fr1", Read(2, 10, "15M"))[0],
			};
		assertEquals("r3", ffactory.findFirstAssociatedSAMRecord(orderedRealigned[0]).getReadName());
		assertEquals("r2", ffactory.findFirstAssociatedSAMRecord(orderedRealigned[1]).getReadName());
		assertEquals(null, ffactory.findFirstAssociatedSAMRecord(orderedRealigned[2]));
		assertEquals("r1", ffactory.findFirstAssociatedSAMRecord(orderedRealigned[3]).getReadName());
		
		assertEquals(null, bfactory.findFirstAssociatedSAMRecord(orderedRealigned[0]));
		assertEquals(null, bfactory.findFirstAssociatedSAMRecord(orderedRealigned[1]));
		assertEquals("r2", bfactory.findFirstAssociatedSAMRecord(orderedRealigned[2]).getReadName());
		assertEquals(null, bfactory.findFirstAssociatedSAMRecord(orderedRealigned[3]));
	}
}
