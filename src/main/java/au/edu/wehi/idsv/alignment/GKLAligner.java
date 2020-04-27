package au.edu.wehi.idsv.alignment;

import com.intel.gkl.smithwaterman.IntelSmithWaterman;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

public class GKLAligner implements Aligner {
	private final IntelSmithWaterman isw;
	private final SWParameters parameters;

	public GKLAligner(int match, int mismatch, int gapOpen, int gapExtend, IntelSmithWaterman isw) {
		this.parameters = new SWParameters(match, mismatch, -gapOpen, -gapExtend);
		this.isw = isw;
	}

	@Override
	public Alignment align_smith_waterman(byte[] seq, byte[] ref) {
		SWNativeAlignerResult result = isw.align(ref, seq, parameters, SWOverhangStrategy.SOFTCLIP);
		Alignment alignment = new Alignment(result.alignment_offset, result.cigar);
		return alignment;
	}
}
