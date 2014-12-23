package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

public class MockDirectedBreakpoint extends DiscordantReadPair {
	private BreakpointSummary breakend;
	private String id;
	public MockDirectedBreakpoint(BreakpointSummary bs) {
		this(bs, null);
	}
	public MockDirectedBreakpoint(BreakpointSummary bs, String id) {
		super(DP(bs)[0], DP(bs)[1], TestHelper.SES(), TestHelper.getContext());
		this.breakend = bs;
		if (id != null) {
			this.id = id;
		}
	}
	private static SAMRecord[] DP(BreakpointSummary bs) {
		return TestHelper.DP(
				bs.referenceIndex,
				bs.start,
				"1M",
				bs.direction == BreakendDirection.Forward,
				bs.referenceIndex2,
				bs.start2,
				"1M",
				bs.direction2 == BreakendDirection.Backward);
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
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
