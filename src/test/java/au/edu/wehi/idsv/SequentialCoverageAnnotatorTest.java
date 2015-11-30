package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;



public class SequentialCoverageAnnotatorTest extends TestHelper {
	public VariantContextDirectedEvidence go(List<SAMRecord> ref, VariantContextDirectedEvidence toAnnotate) {
		Collections.sort(ref, new SAMRecordCoordinateComparator());
		return new SequentialCoverageAnnotator(
				getContext(),
				ImmutableList.of(toAnnotate).iterator(),
				Lists.<ReferenceCoverageLookup>newArrayList(new SequentialReferenceCoverageLookup(ref.iterator(), IDSV(ref), new SAMFlagReadPairConcordanceCalculator(IDSV(ref)), 1024)))
			.annotate(toAnnotate);
	}
	@Test
	public void should_add_reference_counts() {
		VariantContextDirectedEvidence result = go(L(
				RP(0, 1, 100, 5),
				RP(0, 2, 100, 5),
				RP(0, 9, 100, 10)),
			(VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 10, 10, 1, BWD, 100, 100), "")
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
					.breakend(new BreakendSummary(0, FWD, i, i), "")
					.make());
			assertEquals(0, result.getReferenceReadCount(0));
			
			result = go(L(
					Read(0, 2, 1)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, BWD, i, i), "")
					.make());
			assertEquals(0, result.getReferenceReadCount(0));
			
			result = go(L(
					Read(0, 2, 2)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, FWD, i, i), "")
					.make());
			assertEquals(i == 2 ? 1 : 0, result.getReferenceReadCount(0));
			
			result = go(L(
					Read(0, 2, 2)),
				(VariantContextDirectedEvidence)minimalBreakend()
					.breakend(new BreakendSummary(0, BWD, i, i), "")
					.make());
			assertEquals(Integer.toString(i), i == 3 ? 1 : 0, result.getReferenceReadCount(0));
		}
	}
	@Test
	public void should_not_apply_sv_filters() {
		VariantContextDirectedEvidence result = go(new ArrayList<SAMRecord>(),
			(VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 10, 10, 1, BWD, 100, 100), "")
				.make());
		assertFalse(result.isFiltered());
	}
}
