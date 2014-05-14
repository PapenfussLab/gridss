package au.edu.wehi.socrates;

import htsjdk.samtools.SAMRecord;

import com.google.common.collect.PeekingIterator;

public class SequentialBreakendAnnotator {
	private final PeekingIterator<SAMRecord> allReadIt;
	private final PeekingIterator<DirectedBreakpoint> evidenceIt;
	private final ProcessingContext context;
	private final int referenceIndex;
	public SequentialBreakendAnnotator(
			ProcessingContext context,
			int referenceIndex,
			PeekingIterator<SAMRecord> allReads,
			PeekingIterator<DirectedBreakpoint> evidence) {
		this.context = context;
		this.allReadIt = allReads;
		this.evidenceIt = evidence;
		this.referenceIndex = referenceIndex;
	}
	public VariantContextDirectedBreakpoint annotate(VariantContextDirectedBreakpoint variant) {
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(context, variant);
	}
}
