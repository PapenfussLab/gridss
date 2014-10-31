package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

public class SAMEvidenceSourceTest extends IntermediateFilesTest {
	@Test
	public void per_chr_iterator_should_return_all_evidence() {
		createInput(new SAMRecord[] { Read(1, 1, "50M50S") },
				RP(0, 200, 100), // max frag size
				DP(1, 1, "100M", true, 2, 5, "100M", true),
			   DP(1, 2, "100M", true, 2, 4, "100M", true),
			   DP(1, 3, "100M", true, 2, 6, "100M", true),
			   OEA(1, 4, "100M", false));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(true), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(8, list.size()); // 1 SC + 3 * 2 DP + 1 OEA
	}
	@Test
	public void should_stop_metric_calculation_after_max_records() {
		ProcessingContext pc = getCommandlineContext(true);
		pc.setCalculateMetricsRecordCount(2);
		createInput(RP(0, 100, 200, 100), RP(0, 400, 600, 100));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, false);
		source.completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS));
		assertEquals(200-100+100, source.getMetrics().getMaxFragmentSize());
		pc.setCalculateMetricsRecordCount(1000);
		
		source = new SAMEvidenceSource(pc, input, false);
		source.completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS));
		assertEquals(600-400+100, source.getMetrics().getMaxFragmentSize());
	}
	@Test
	public void per_chr_iterator_should_iterator_over_chr_in_dictionary_order() {
		createInput(
				Read(0, 1, "50M50S"),
				Read(1, 1, "50M50S"),
				Read(2, 1, "50M50S"),
				Read(3, 1, "50M50S"),
				Read(4, 1, "50M50S")
				);
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(true), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		for (int i = 0; i <= 4; i++) {
			assertEquals(i, list.get(i).getBreakendSummary().referenceIndex);
		}
	}
	@Test
	public void per_chr_iterator_chr_should_return_only_evidence() {
		createInput(
				Read(0, 1, "50M50S"),
				Read(1, 1, "50M50S"),
				Read(2, 1, "50M50S"),
				Read(3, 1, "50M50S"),
				Read(4, 1, "50M50S")
				);
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(true), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator("polyACGT"));
		assertEquals(1, list.size());
		assertEquals(1, list.get(0).getBreakendSummary().referenceIndex);
	}
	@Test
	public void iterator_should_return_all_evidence() {
		createInput(new SAMRecord[] { Read(1, 1, "50M50S") },
				RP(0, 100, 200, 100),
				DP(1, 1, "100M", true, 2, 5, "100M", true),
			   DP(1, 2, "100M", true, 2, 4, "100M", true),
			   DP(1, 3, "100M", true, 2, 6, "100M", true),
			   OEA(1, 4, "100M", false));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(8, list.size()); // 1 SC + 3 * 2 DP + 1 OEA
	}
	@Test
	public void iterator_should_iterator_over_chr_in_dictionary_order() {
		createInput(
				Read(0, 1, "50M50S"),
				Read(1, 1, "50M50S"),
				Read(2, 1, "50M50S"),
				Read(3, 1, "50M50S"),
				Read(4, 1, "50M50S")
				);
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		for (int i = 0; i <= 4; i++) {
			assertEquals(i, list.get(i).getBreakendSummary().referenceIndex);
		}
	}
	@Test
	public void iterator_chr_should_return_only_evidence() {
		createInput(
				Read(0, 1, "50M50S"),
				Read(1, 1, "50M50S"),
				Read(2, 1, "50M50S"),
				Read(3, 1, "50M50S"),
				Read(4, 1, "50M50S")
				);
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator("polyACGT"));
		assertEquals(1, list.size());
		assertEquals(1, list.get(0).getBreakendSummary().referenceIndex);
	}
	@Test
	public void should_set_evidence_source_to_self() {
		createInput(Read(0, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(1, list.size());
		assertEquals(source, list.get(0).getEvidenceSource());
	}
	@Test
	public void should_default_fragment_size_to_read_length_for_unpaired_reads() {
		createInput(Read(0, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		assertEquals(100, source.getMetrics().getMaxFragmentSize());
	}
	@Test
	public void iterator_evidence_should_be_sorted_by_evidence_natural_ordering() {
		createInput(new SAMRecord[] { Read(1, 1, "50M50S") },
				new SAMRecord[] { Read(1, 1, "50S50M") },
				RP(0, 100, 200, 100),
				DP(1, 1, "100M", true, 2, 5, "100M", true),
			   DP(1, 2, "100M", false, 2, 4, "100M", true),
			   DP(1, 3, "100M", true, 2, 6, "100M", true),
			   OEA(1, 4, "100M", false),
			   OEA(1, 4, "100M", true));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		List<DirectedEvidence> sorted = Lists.newArrayList(result);
		Collections.sort(sorted, DirectedEvidenceOrder.ByNatural);
		assertEquals(sorted, result);
	}
	@Test
	public void iterator_should_match_soft_clip_realignment_with_soft_clip() {
		createInput(
				withReadName("r1", Read(0, 1, "15M15S")),
				withReadName("r2", Read(1, 2, "15M15S")),
				withReadName("r3", Read(2, 3, "15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(getCommandlineContext().getFileSystemContext().getRealignmentBam(input), SortOrder.unsorted, 
				withReadName("0#1#fr1", Read(2, 10, "15M"))[0],
				withReadName("1#2#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#fr3", Read(0, 10, "15M"))[0]
		);
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(3, result.size());
		assertTrue(result.get(0) instanceof RealignedSoftClipEvidence);
		assertTrue(result.get(1) instanceof RealignedSoftClipEvidence);
		assertTrue(result.get(2) instanceof RealignedSoftClipEvidence);
	}
	@Test
	public void iterator_should_iterator_over_both_forward_and_backward_soft_clips() {
		createInput(withReadName("r1", Read(0, 1, "15S15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(2, result.size());
	}
	@Test
	public void iterator_should_realign_both_forward_and_backward_soft_clips() {
		createInput(
				withReadName("r1", Read(0, 1, "15S15M15S")),
				withReadName("r2", Read(1, 2, "15S15M15S")),
				withReadName("r3", Read(2, 3, "15S15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(getCommandlineContext().getFileSystemContext().getRealignmentBam(input), SortOrder.unsorted, 
				withReadName("0#1#fr1", Read(2, 10, "15M"))[0],
				withReadName("0#1#br1", Read(1, 15, "15M"))[0],
				withReadName("1#2#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#fr3", Read(0, 10, "15M"))[0],
				withReadName("2#3#br3", Read(0, 100, "15M"))[0]
		);
		List<RealignedSoftClipEvidence> result = Lists.newArrayList(Iterators.filter(source.iterator(), RealignedSoftClipEvidence.class));
		
		assertEquals(5, result.size());
	}
	@Test
	public void iterator_should_include_remote_soft_clips() {
		createInput(
				withReadName("r1", Read(0, 1, "15S15M15S")),
				withReadName("r2", Read(1, 2, "15S15M15S")),
				withReadName("r3", Read(2, 3, "15S15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, false);
		source.completeSteps(ProcessStep.ALL_STEPS);
		createBAM(getCommandlineContext().getFileSystemContext().getRealignmentBam(input), SortOrder.unsorted, 
				withReadName("0#1#fr1", Read(2, 10, "15M"))[0],
				withReadName("0#1#br1", Read(1, 15, "15M"))[0],
				withReadName("1#2#fr2", Read(1, 10, "15M"))[0],
				withReadName("2#3#fr3", Read(0, 10, "15M"))[0],
				withReadName("2#3#br3", Read(0, 100, "15M"))[0]);
		source.completeSteps(ProcessStep.ALL_STEPS);
		
		List<RealignedRemoteSoftClipEvidence> remote = Lists.newArrayList(Iterators.filter(source.iterator(), RealignedRemoteSoftClipEvidence.class));
		assertEquals(5, remote.size());
		
		List<RealignedSoftClipEvidence> result = Lists.newArrayList(Iterators.filter(source.iterator(), RealignedSoftClipEvidence.class));
		assertEquals(5, result.size() - remote.size());
	}
}
