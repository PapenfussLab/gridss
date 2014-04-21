package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

public class DeBruijnAssemblyReadProcessor implements SequentialDirectedBreakpointEvidenceProcessor {

	@Override
	public Iterable<DirectedBreakpoint> getEvidenceAtPosition(int position) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void addSoftclipSupportingForwardBreakpoint(SAMRecord record) {
		// TODO Auto-generated method stub

	}

	@Override
	public void addSoftclipSupportingBackwardBreakpoint(SAMRecord record) {
		// TODO Auto-generated method stub

	}

	@Override
	public void addReadPairSupportingForwardBreakpoint(NonReferenceReadPair pair) {
		// TODO Auto-generated method stub

	}

	@Override
	public void addReadPairSupportingBackwardBreakpoint(
			NonReferenceReadPair pair) {
		// TODO Auto-generated method stub

	}

}
