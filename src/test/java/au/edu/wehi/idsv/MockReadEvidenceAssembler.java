package au.edu.wehi.idsv;

import java.util.List;

import com.google.common.collect.ImmutableList;

public class MockReadEvidenceAssembler implements ReadEvidenceAssembler {
	private List<SAMRecordAssemblyEvidence> evidence;

	public MockReadEvidenceAssembler(List<SAMRecordAssemblyEvidence> evidence) {
		this.evidence = evidence;
	}

	@Override
	public Iterable<SAMRecordAssemblyEvidence> addEvidence(DirectedEvidence evidence) {
		return ImmutableList.of();
	}

	@Override
	public Iterable<SAMRecordAssemblyEvidence> endOfEvidence() {
		return evidence;
	}

	@Override
	public String getStateSummaryMetrics() {
		return "mock";
	}
}
