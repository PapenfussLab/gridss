package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class AggregateReferenceCoverageLookup implements ReferenceCoverageLookup {
	private List<ReferenceCoverageLookup> underlying;
	public AggregateReferenceCoverageLookup(Collection<ReferenceCoverageLookup> underlying) {
		this.underlying = new ArrayList<ReferenceCoverageLookup>(underlying);
	}
	@Override
	public int readsSupportingNoBreakendAfter(int referenceIndex, int position) {
		int count = 0;
		for (ReferenceCoverageLookup lookup : underlying) {
			count += lookup.readsSupportingNoBreakendAfter(referenceIndex, position);
		}
		return count;
	}

	@Override
	public int readPairsSupportingNoBreakendAfter(int referenceIndex, int position) {
		int count = 0;
		for (ReferenceCoverageLookup lookup : underlying) {
			count += lookup.readPairsSupportingNoBreakendAfter(referenceIndex, position);
		}
		return count;
	}

}
