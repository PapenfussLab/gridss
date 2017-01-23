package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;



public class SequentialCoverageAnnotatorTest extends TestHelper {
	@SuppressWarnings("resource")
	public VariantContextDirectedEvidence go(List<SAMRecord> ref, VariantContextDirectedEvidence toAnnotate) {
		Collections.sort(ref, new SAMRecordCoordinateComparator());
		return new SequentialCoverageAnnotator<VariantContextDirectedEvidence>(
				getContext(),
				ImmutableList.of(toAnnotate).iterator(),
				Lists.<ReferenceCoverageLookup>newArrayList(new SequentialReferenceCoverageLookup(ref.iterator(), IDSV(ref), new SAMFlagReadPairConcordanceCalculator(IDSV(ref)), 1024, 0)),
				MoreExecutors.newDirectExecutorService())
			.annotate(toAnnotate);
	}
	@Test
	public void should_add_reference_counts() {
		VariantContextDirectedEvidence result = go(L(
				RP(0, 1, 100, 5),
				RP(0, 2, 100, 5),
				RP(0, 9, 100, 10)),
			(VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 10, 1, BWD, 100), "")
				.make());
		assertEquals(1, result.getReferenceReadCount(0));
		assertEquals(2, result.getReferenceReadPairCount(0));
	}
	@Test
	public void should_count_reference_support_spanning_breakend_position() {
		VariantContextDirectedEvidence result;
		for (int i = 1; i < 10; i++) {
			result = go(L(
					Read(0, 2, 1)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, FWD, i), "")
					.make());
			assertEquals(0, result.getReferenceReadCount(0));
			
			result = go(L(
					Read(0, 2, 1)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, BWD, i), "")
					.make());
			assertEquals(0, result.getReferenceReadCount(0));
			
			result = go(L(
					Read(0, 2, 2)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, FWD, i), "")
					.make());
			assertEquals(i == 2 ? 1 : 0, result.getReferenceReadCount(0));
			
			result = go(L(
					Read(0, 2, 2)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, BWD, i), "")
					.make());
			assertEquals(Integer.toString(i), i == 3 ? 1 : 0, result.getReferenceReadCount(0));
		}
	}
	@Test
	public void should_not_apply_sv_filters() {
		VariantContextDirectedEvidence result = go(new ArrayList<SAMRecord>(),
			(VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 10, 1, BWD, 100), "")
				.make());
		assertFalse(result.isFiltered());
	}
	@Test
	public void should_use_homology_interval_minimum_as_reference_coverage() {
		//  1234567890 >
		//  MMMMMMMMMMM
		//  MMMMM
		//    MMMMMMMMMM
		//	  MMMMMMMMMMM
		//	MM
		VariantContextDirectedEvidence result = go(new ArrayList<>(ImmutableList.of(
				Read(0, 1, "11M"),
				Read(0, 1, "5M"),
				Read(0, 3, "10M"),
				Read(0, 3, "11M"),
				Read(0, 1, "2M")
			)), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 1, 1, 10, 1, BWD, 100, 100, 100), "")
				.make());
		assertEquals(2, result.getReferenceReadCount());
	}
	@Test
	public void category_count_should_match_ProcessContext_count() {
		VariantContextDirectedEvidence result = go(new ArrayList<>(ImmutableList.of(
				Read(0, 1, "11M"),
				Read(0, 1, "5M"),
				Read(0, 3, "10M"),
				Read(0, 3, "11M"),
				Read(0, 1, "2M")
			)), (VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 1, 1, 10, 1, BWD, 100, 100, 100), "")
				.make());
		assertEquals(2, ((int[])result.getAttribute("REF")).length);
	}
}
