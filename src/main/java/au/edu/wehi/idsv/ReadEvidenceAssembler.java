package au.edu.wehi.idsv;

public interface ReadEvidenceAssembler {
	Iterable<VariantContextDirectedEvidence> addEvidence(DirectedEvidence evidence);
	Iterable<VariantContextDirectedEvidence> endOfEvidence();
	String getStateSummaryMetrics();
}
