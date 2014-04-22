package au.edu.wehi.socrates;

public interface ReadEvidenceAssembler {
	Iterable<AssemblyEvidence> addEvidence(DirectedEvidence evidence);
	Iterable<AssemblyEvidence> endOfEvidence();
}
