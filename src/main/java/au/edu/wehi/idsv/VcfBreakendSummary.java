package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfConstants;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.ImmutableList;

/**
 * Extracts breakend/breakpoint information from a VCF variant
 * @author cameron.d
 *
 */
public class VcfBreakendSummary {
	public final BreakendSummary location;
	public final String breakpointSequence;
	public final String anchorSequence;
	public VcfBreakendSummary(ProcessingContext processContext, VariantContext variant) {
		List<Allele> altList = variant.getAlternateAlleles();
		if (altList.size() != 1) {
			location = null;
			breakpointSequence = null;
			anchorSequence = null;
			return;
		}
		String alt = variant.getAlternateAllele(0).getDisplayString();
		alt = alt.replace(VcfConstants.VCF41BREAKEND_REPLACEMENT, VcfConstants.VCF42BREAKEND);
		if (variant.getReference().length() >= alt.length()) {
			// Technically this is valid (eg: {"AAA", "A."} = breakend with deletion), we just can't handle these yet
			location = null;
			breakpointSequence = null;
			anchorSequence = null;
			return;
		}
		if (alt.length() < 2) {
			location = null;
			breakpointSequence = null;
			anchorSequence = null;
			return;
		}
		BreakendDirection direction, remoteDirection = null;
		String localSequence;
		String remoteContig = null;
		if (alt.charAt(0) == '.') {
			// .BreakpointReference
			direction = BreakendDirection.Backward;
			localSequence = alt.substring(1);
		} else if (alt.charAt(alt.length() - 1) == '.') {
			// ReferenceBreakpoint.
			direction = BreakendDirection.Forward;
			localSequence = alt.substring(0, alt.length() - 1);
		} else if (alt.charAt(0) == '[' || alt.charAt(0) == ']') {
			// [Remote[BreakpointReference
			direction = BreakendDirection.Backward;
			remoteDirection = alt.charAt(0) == ']' ? BreakendDirection.Forward : BreakendDirection.Backward;
			String[] split = alt.split("[\\[\\]]");
			remoteContig = split[1];
			localSequence = split[2];
		} else if (alt.charAt(alt.length() - 1) == '[' || alt.charAt(alt.length() - 1) == ']') {
			// ReferenceBreakpoint[Remote[
			direction = BreakendDirection.Forward;
			remoteDirection = alt.charAt(alt.length() - 1) == ']' ? BreakendDirection.Forward : BreakendDirection.Backward;
			String[] split = alt.split("[\\[\\]]");
			remoteContig = split[1];
			localSequence = split[0];
		} else {
			// not breakend!
			location = null;
			breakpointSequence = null;
			anchorSequence = null;
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
		int refLength = variant.getReference().length();
		int localPosition;
		if (direction == BreakendDirection.Forward) {
			localPosition = variant.getEnd();
			// anchor - breakpoint
			anchorSequence = localSequence.substring(0, refLength);
			breakpointSequence = localSequence.substring(anchorSequence.length());
		} else {
			localPosition = variant.getStart();
			// breakpoint - anchor
			breakpointSequence = localSequence.substring(0, localSequence.length() - refLength);
			anchorSequence = localSequence.substring(breakpointSequence.length());
		}
		int ciStart = 0, ciEnd = 0;
		int rciStart = 0, rciEnd = 0;
		if (variant.hasAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)) {
			List<Integer> ci = AttributeConverter.asIntList(variant.getAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)); 
			if (ci.size() == 2) {
				ciStart = ci.get(0);
				ciEnd = ci.get(1);
			} else {
				throw new IllegalStateException(String.format("Error parsing attribute %s of %s. Expected 2 integer values, found %d", VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, super.toString(), ci.size()));
			}
		}
		if (variant.hasAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())) {
			List<Integer> rci = AttributeConverter.asIntList(variant.getAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())); 
			if (rci.size() == 2) {
				rciStart = rci.get(0);
				rciEnd = rci.get(1);
			} else {
				throw new IllegalStateException(String.format("Error parsing attribute %s of %s. Expected 2 integer values, found %d", VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute(), super.toString(), rci.size()));
			}
		}
		if (remoteDirection != null) {
			location = new BreakpointSummary(IdsvVariantContext.getReferenceIndex(processContext, variant), direction, localPosition + ciStart, localPosition + ciEnd,
					processContext.getDictionary().getSequenceIndex(remoteContig), remoteDirection, remotePosition + rciStart, remotePosition + rciEnd);
		} else {
			location = new BreakendSummary(IdsvVariantContext.getReferenceIndex(processContext, variant), direction, localPosition + ciStart, localPosition + ciEnd);
		}
	}
}
