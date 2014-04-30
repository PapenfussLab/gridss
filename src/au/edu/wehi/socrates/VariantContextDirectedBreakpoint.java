package au.edu.wehi.socrates;

import java.nio.charset.StandardCharsets;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;

/**
 * VCF Breakend record
 * see Section 5.4.9 of http://samtools.github.io/hts-specs/VCFv4.2.pdf for details of breakends
 * @author Daniel Cameron
 *
 */
public class VariantContextDirectedBreakpoint extends SocratesVariantContext implements DirectedBreakpoint {
	private BreakpointLocation location;
	private String breakpointSequence;
	private String anchorSequence;
	protected VariantContextDirectedBreakpoint(ProcessingContext processContext, VariantContext context) {
		super(processContext, context);
		recalcFields();
	}
	private void recalcFields() {
		location = null;
		breakpointSequence = null;
		anchorSequence = null;
		
		List<Allele> altList = getAlternateAlleles();
		if (altList.size() != 1) return;
		String alt = getAlternateAllele(0).getDisplayString();
		if (getReference().length() >= alt.length()) {
			// Technically this is valid (eg: {"AAA", "A."} = breakend with deletion), we just can't handle these yet
			return;
		}
		if (alt.length() < 2) return;
		BreakpointDirection direction, remoteDirection = null;
		String localSequence;
		String remoteContig = null;
		if (alt.charAt(0) == '.') {
			// .BreakpointReference
			direction = BreakpointDirection.Backward;
			localSequence = alt.substring(1);
		} else if (alt.charAt(alt.length() - 1) == '.') {
			// ReferenceBreakpoint.
			direction = BreakpointDirection.Forward;
			localSequence = alt.substring(0, alt.length() - 1);
		} else if (alt.charAt(0) == '[' || alt.charAt(0) == ']') {
			// [Remote[BreakpointReference
			direction = BreakpointDirection.Backward;
			remoteDirection = alt.charAt(0) == '[' ? BreakpointDirection.Forward : BreakpointDirection.Backward;
			String[] split = alt.split("[\\[\\]]");
			remoteContig = split[1];
			localSequence = split[2];
		} else if (alt.charAt(alt.length() - 1) == '[' || alt.charAt(alt.length() - 1) == ']') {
			// ReferenceBreakpoint[Remote[
			direction = BreakpointDirection.Forward;
			remoteDirection = alt.charAt(alt.length() - 1) == '[' ? BreakpointDirection.Forward : BreakpointDirection.Backward;
			String[] split = alt.split("[\\[\\]]");
			remoteContig = split[1];
			localSequence = split[0];
		} else {
			// not breakend!
			return;
		}
		int remotePosition = 0;
		if (StringUtils.isNotEmpty(remoteContig)) {
			// flanking square brackets have already been removed
			// format of chr:pos so breakend should always specify a contig position
			String[] components = remoteContig.split(":");
			remoteContig = components[0];
			if (components.length > 1) {
				remotePosition = Integer.parseInt(components[1]);
			}
		}
		int refLength = getReference().length();
		int localPosition;
		if (direction == BreakpointDirection.Forward) {
			localPosition = getEnd();
			// anchor - breakpoint
			anchorSequence = localSequence.substring(0, refLength);
			breakpointSequence = localSequence.substring(anchorSequence.length());
		} else {
			localPosition = getStart();
			// breakpoint - anchor
			breakpointSequence = localSequence.substring(0, localSequence.length() - refLength);
			anchorSequence = localSequence.substring(breakpointSequence.length());
		}
		if (remoteDirection != null) {
			location = new BreakpointInterval(getReferenceIndex(), direction, localPosition, localPosition,
					processContext.getDictionary().getSequenceIndex(remoteContig), remoteDirection, remotePosition, remotePosition,
					getPhredScaledQual());
		} else {
			location = new BreakpointLocation(getReferenceIndex(), direction, localPosition, localPosition,
					getPhredScaledQual());
		}
	}
	@Override
	public BreakpointLocation getBreakpointLocation() {
		if (location == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return location;
	}
	@Override
	public String getEvidenceID() {
		return getID();
	}
	@Override
	public byte[] getBreakpointSequence() {
		return breakpointSequence.getBytes(StandardCharsets.US_ASCII);
	}
	public String getBreakpointSequenceString() {
		if (breakpointSequence == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakpointSequence;
	}
	public String getAnchorSequenceString() {
		if (anchorSequence == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return anchorSequence;
	}
	@Override
	public byte[] getBreakpointQuality() {
		return null;
	}
	@Override
	public boolean isValid() {
		return location != null;
	}
}
