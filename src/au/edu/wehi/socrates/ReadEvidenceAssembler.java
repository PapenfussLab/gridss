package au.edu.wehi.socrates;

public interface ReadEvidenceAssembler {
	Iterable<DirectedBreakpointAssembly> addEvidence(DirectedEvidence evidence);
	Iterable<DirectedBreakpointAssembly> endOfEvidence();
}
