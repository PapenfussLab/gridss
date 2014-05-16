package au.edu.wehi.socrates;

public interface ReadEvidenceAssembler {
	Iterable<VariantContextDirectedBreakpoint> addEvidence(DirectedEvidence evidence);
	Iterable<VariantContextDirectedBreakpoint> endOfEvidence();
}
