package au.edu.wehi.idsv;

public interface ReadEvidenceAssembler {
	Iterable<AssemblyEvidence> addEvidence(DirectedEvidence evidence);
	Iterable<AssemblyEvidence> endOfEvidence();
	String getStateSummaryMetrics();
}
