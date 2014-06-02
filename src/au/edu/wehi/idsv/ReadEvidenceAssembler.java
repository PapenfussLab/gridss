package au.edu.wehi.idsv;

public interface ReadEvidenceAssembler {
	Iterable<VariantContextDirectedBreakpoint> addEvidence(DirectedEvidence evidence);
	Iterable<VariantContextDirectedBreakpoint> endOfEvidence();
}
