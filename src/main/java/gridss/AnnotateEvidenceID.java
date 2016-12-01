package gridss;

import java.util.Iterator;
import java.util.stream.Collectors;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.SequentialEvidenceAllocator;
import au.edu.wehi.idsv.SequentialEvidenceAllocator.VariantEvidenceSupport;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import picard.cmdline.CommandLineProgramProperties;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
        usage = "Annotates variant calls with the GRIDSS evidenceIDs of the supporting evidence.",  
        usageShort = "Annotates variant calls with the GRIDSS evidenceIDs of the supporting evidence."
)
public class AnnotateEvidenceID extends VcfTransformCommandLineProgram {
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls) {
		CloseableIterator<DirectedEvidence> evidence = getEvidenceIterator();
		Iterator<VariantEvidenceSupport> annotator = new SequentialEvidenceAllocator(getContext(), calls, evidence, maximumWindowSize(), false);
		Iterator<VariantContextDirectedBreakpoint> it = Iterators.transform(annotator, ves -> annotate(ves));
		return new AutoClosingIterator<>(it, evidence, calls);
	}
	private VariantContextDirectedBreakpoint annotate(VariantEvidenceSupport ves) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), ves.variant);
		builder.attribute(VcfAttributes.EVIDENCE_ID, ves.support.stream().map(e -> e.getEvidenceID()).collect(Collectors.toList()));
		IdsvVariantContext variant = builder.make();
		return (VariantContextDirectedBreakpoint)variant;
	}
}
