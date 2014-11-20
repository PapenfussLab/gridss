package au.edu.wehi.idsv;

import java.util.List;

import com.google.common.collect.ImmutableList;

public class MockReadEvidenceAssembler implements ReadEvidenceAssembler {
	private List<AssemblyEvidence> evidence;

	public MockReadEvidenceAssembler(List<AssemblyEvidence> evidence) {
		this.evidence = evidence;
	}

	@Override
	public Iterable<AssemblyEvidence> addEvidence(DirectedEvidence evidence) {
		return ImmutableList.of();
	}

	@Override
	public Iterable<AssemblyEvidence> endOfEvidence() {
		return evidence;
	}

	@Override
	public String getStateSummaryMetrics() {
		return "mock";
	}
}
