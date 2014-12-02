package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;



public class SequentialEvidenceAnnotatorTest extends TestHelper {
	public VariantContextDirectedEvidence go(List<DirectedEvidence> evidence, VariantContextDirectedEvidence toAnnotate) {
		return new SequentialEvidenceAnnotator(getContext(), ImmutableList.of(toAnnotate).iterator(), Iterators.peekingIterator(evidence.iterator())).annotate(toAnnotate);
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
		SAMRecordAssemblyEvidence be = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 1, B("TGGAAT"), B("TAAAAT"), 0, 0);
		// Untemplated sequence
		SAMRecord r = Read(0, 10, "2S3M");
		r.setReadBases(B("GGAAT"));
		RealignedSAMRecordAssemblyEvidence ass = (RealignedSAMRecordAssemblyEvidence) AssemblyFactory.incorporateRealignment(getContext(), be, r);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 0, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)ass
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals("AGG[polyA:10[", result.getAlternateAllele(0).getDisplayString());
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
}
