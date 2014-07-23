package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.Collections;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes.Subset;

import com.google.common.collect.Iterators;



public class SequentialBreakendAnnotatorTest extends TestHelper {
	public VariantContextDirectedEvidence go(List<DirectedEvidence> evidence, VariantContextDirectedEvidence toAnnotate) {
		return new SequentialBreakendAnnotator(getContext(), null, null, Iterators.peekingIterator(evidence.iterator())).annotate(toAnnotate);
	}
	public VariantContextDirectedEvidence go(List<DirectedEvidence> evidence, List<SAMRecord> ref, VariantContextDirectedEvidence toAnnotate) {
		Collections.sort(ref, new SAMRecordCoordinateComparator());
		return new SequentialBreakendAnnotator(getContext(), new SequentialReferenceCoverageLookup(ref.iterator(), 1000), null, Iterators.peekingIterator(evidence.iterator())).annotate(toAnnotate);
	}
	@Test
	public void should_add_breakpoint_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(1, result.getEvidenceCountSoftClip(null));
	}
	@Test
	public void should_add_breakend_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(1, result.getEvidenceCountSoftClip(null));
	}
	@Test
	public void should_use_breakpoint_assembly_sequence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		IdsvVariantContextBuilder assBuilder = new IdsvVariantContextBuilder(getContext(), AE());
		assBuilder.breakpoint(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 10, 10), "TT");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)assBuilder.make()
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals("ATT[polyA:10[", result.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_merge_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M")),
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(2, result.getEvidenceCountSoftClip(null));
	}
	@Test
	public void should_not_add_breakend_nonsupporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), BWD, Read(0, 1, "3S1M"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(0, result.getEvidenceCountSoftClip(null));
	}
	@Test
	public void should_not_add_breakpoint_nonsupporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("TTTTTTTT", Read(0, 1, "3S1M3S"))[0], withSequence("TTT", Read(0, 12, "3M"))[0])
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(0, result.getEvidenceCountSoftClip(null));
	}
	@Test
	public void should_add_reference_counts() {
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("TTTTTTTT", Read(0, 1, "3S1M3S"))[0], withSequence("TTT", Read(0, 12, "3M"))[0])),
			L(
				RP(0, 1, 100, 5),
				RP(0, 2, 100, 5),
				RP(0, 9, 100, 10)),
			(VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 10, 10, 1, BWD, 100, 100), "")
				.make());
		assertEquals(1, result.getReferenceReadCount(Subset.NORMAL));
		assertEquals(2, result.getReferenceReadPairCount(Subset.NORMAL));
	}
}
