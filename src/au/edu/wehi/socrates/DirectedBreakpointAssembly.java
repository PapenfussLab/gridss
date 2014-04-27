package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

public class DirectedBreakpointAssembly extends VariantContextDirectedBreakpoint {
	public static final String SOURCE_NAME = "socrates";
	public static DirectedBreakpointAssembly create(DirectedBreakpointAssembly variant, SAMRecord realignment) {
		throw new IllegalStateException("not yet implemented");
	}
	public static DirectedBreakpointAssembly create(
			SAMSequenceDictionary dictionary,
			ReferenceSequenceFile reference,
			String assemblerName,
			int referenceIndex,
			int position,
			BreakpointDirection direction,
			byte[] breakpointSequence,
			byte[] fullAssembly,
			double breakpointQuality
			) {
		String chr = dictionary.getSequence(referenceIndex).getSequenceName();
		VariantContextBuilder builder = new VariantContextBuilder()
			.source(SOURCE_NAME)
			.id(String.format("%s-%s:%d-%s", assemblerName, chr, position, direction == BreakpointDirection.Forward ? "f" : "b"))
			.chr(chr)
			.start(position)
			.stop(position)
			.log10PError(breakpointQuality)
			.attributes(null)
			.attribute(VcfConstants.ASSEMBLY_PROGRAM, assemblerName)
			.attribute(VcfConstants.ASSEMBLY_CONSENSUS, new String(fullAssembly, StandardCharsets.US_ASCII));
		String referenceBase = new String(reference.getSubsequenceAt(chr, position, position).getBases(), StandardCharsets.US_ASCII);
		String alt;
		String breakStr = new String(breakpointSequence, StandardCharsets.US_ASCII);
		// TODO: use reference bases or anchor assembly?
		if (direction == BreakpointDirection.Forward) {
			alt = referenceBase + breakStr + '.';
		} else {
			alt = '.' + breakStr + referenceBase;
		}
		builder.alleles(referenceBase, alt);
		return new DirectedBreakpointAssembly(dictionary, builder.make());
	}
	protected DirectedBreakpointAssembly(SAMSequenceDictionary dictionary, VariantContext context) {
		super(dictionary, context);
	}
	@Override
	public boolean isValid() {
		return super.isValid() && getAssemblerProgram() != null;
	};
	public String getAssemblerProgram() { return getAttributeAsString(VcfConstants.ASSEMBLY_PROGRAM, null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfConstants.ASSEMBLY_CONSENSUS, ""); }
	public String getAnchorConsensus() {
		String consensus = getAssemblyConsensus();
		if (getBreakpointLocation().direction == BreakpointDirection.Forward) {
			return consensus.substring(0, consensus.length() - getBreakpointSequenceString().length());
		} else {
			return consensus.substring(getBreakpointSequenceString().length());
		}
	}
	/**
	 * Number of reads contributing to consensus
	 * @return
	 */
	public int getConsensusReadCount() {
		return getAttributeAsInt(VcfConstants.ASSEMBLY_CONSENSUS_READ_COUNT, 0);
	}
}
