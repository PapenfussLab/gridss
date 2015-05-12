package au.edu.wehi.idsv;

public class MockDirectedEvidence extends UnmappedMateReadPair {
	public BreakendSummary breakend;
	public String id = toString();
	public MockDirectedEvidence(int referenceIndex, BreakendDirection direction, int start) {
		this(new BreakendSummary(referenceIndex, direction, start, start));
	}
	public MockDirectedEvidence(int referenceIndex, BreakendDirection direction, int start, int end) {
		this(new BreakendSummary(referenceIndex, direction, start, end));
	}
	public MockDirectedEvidence(int referenceIndex, int start, String id) {
		this(new BreakendSummary(referenceIndex, BreakendDirection.Forward, start, start), id);
	}
	public MockDirectedEvidence(BreakendSummary bs) {
		this(bs, bs.toString().replace(' ', '_'));
	}
	public MockDirectedEvidence(BreakendSummary bs, String id) {
		super(TestHelper.OEA(0, 1, "1M", true)[0], TestHelper.OEA(0, 1, "1M", true)[1], TestHelper.SES());
		this.breakend = bs;
		if (id != null) {
			this.id = id;
		}
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		return breakend;
	}

	@Override
	public byte[] getBreakendSequence() {
		return null;
	}

	@Override
	public byte[] getBreakendQuality() {
		return null;
	}

	@Override
	public String getEvidenceID() {
		return id;
	}
	@Override
	public SAMEvidenceSource getEvidenceSource() {
		return TestHelper.SES();
	}

	@Override
	public int getLocalMapq() {
		return 1;
	}

	@Override
	public int getLocalBaseLength() {
		return 1;
	}
	@Override
	public int getLocalMaxBaseQual() {
		return 1;
	}

	@Override
	public int getLocalTotalBaseQual() {
		return 1;
	}
}