package au.edu.wehi.socrates;

public interface ReadEvidenceAssembler {
	Iterable<VariantContextEvidence> addEvidence(DirectedEvidence evidence);
	Iterable<VariantContextEvidence> endOfEvidence();
}
