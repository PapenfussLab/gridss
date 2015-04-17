package au.edu.wehi.idsv;

public interface ReadEvidenceAssembler {
	Iterable<SAMRecordAssemblyEvidence> addEvidence(DirectedEvidence evidence);
	Iterable<SAMRecordAssemblyEvidence> endOfEvidence();
	String getStateSummaryMetrics();
}
