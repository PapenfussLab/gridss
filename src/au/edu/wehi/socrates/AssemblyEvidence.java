package au.edu.wehi.socrates;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.variant.variantcontext.VariantContext;

import au.edu.wehi.socrates.util.SAMRecordSummary;

public class AssemblyEvidence extends VariantContextDirectedBreakpoint {
	public static final char BREAKPOINT_ID_CHAR = 'a';
	public AssemblyEvidence(VariantContext context) {
		super(context);
	}
	public boolean isValid() {
		// TODO: check that is this correct
		return true;
	}
}
