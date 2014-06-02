package au.edu.wehi.idsv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfConstants;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;

/**
 * VCF Breakend record
 * see Section 5.4.9 of http://samtools.github.io/hts-specs/VCFv4.2.pdf for details of breakends
 * @author Daniel Cameron
 *
 */
public class VariantContextDirectedBreakpoint extends IdsvVariantContext implements DirectedBreakpoint {
	private final BreakendSummary location;
	private final String breakpointSequence;
	private final String anchorSequence;
	private final byte[] breakpointQual;
	//private static Log LOG = Log.getInstance(VariantContextDirectedBreakpoint.class);
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, VariantContext context) {
		this(processContext, context, null);
	}
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, VariantContext context, byte[] breakpointQual) {
		super(processContext, context);
		this.breakpointQual = breakpointQual;
		// calculate fields
		List<Allele> altList = getAlternateAlleles();
		if (altList.size() != 1) {
			location = null;
			breakpointSequence = null;
			anchorSequence = null;
			return;
		}
		String alt = getAlternateAllele(0).getDisplayString();
		alt = alt.replace(VcfConstants.VCF41BREAKEND_REPLACEMENT, VcfConstants.VCF42BREAKEND);
		if (getReference().length() >= alt.length()) {
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
		int refLength = getReference().length();
		int localPosition;
		if (direction == BreakendDirection.Forward) {
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
		int ciStart = 0, ciEnd = 0;
		if (hasAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)) {
			List<Integer> ci = parseIntList(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY);
			if (ci.size() == 2) {
				ciStart = ci.get(0);
				ciEnd = ci.get(1);
			} else {
				throw new IllegalStateException(String.format("Error parsing attribute %s of %s. Expected 2 integer values, found %d", VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, super.toString(), ci.size()));
			}
		}
		if (remoteDirection != null) {
			location = new BreakpointSummary(getReferenceIndex(), direction, localPosition - ciStart, localPosition + ciEnd,
					processContext.getDictionary().getSequenceIndex(remoteContig), remoteDirection, remotePosition - ciStart, remotePosition + ciEnd,
					getEvidence());
		} else {
			location = new BreakendSummary(getReferenceIndex(), direction, localPosition - ciStart, localPosition + ciEnd,
					getEvidence());
		}
	}
	private EvidenceMetrics getEvidence() {
		EvidenceMetrics m = new EvidenceMetrics();
		for (VcfAttributes a : VcfAttributes.evidenceValues()) {
			m.set(a, getAttributeAsInt(a.attribute(), 0));
		}
		return m;
	}
	@Override
	public BreakendSummary getBreakendSummary() {
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
		return breakpointQual;
	}
	@Override
	public boolean isValid() {
		return location != null;
	}
	public String getAssemblerProgram() { return getAttributeAsString(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), ""); }
	public double getAssemblyQuality() { return getAttributeAsDouble(VcfAttributes.ASSEMBLY_QUALITY.attribute(), 0); }
	/**
	 * Returns an iterator containing only the breakend variants from the given iterator
	 * @param context processing context
	 * @param variants input variants
	 * @return breakend variants
	 */
	public static Iterator<VariantContextDirectedBreakpoint> breakendIterator(final ProcessingContext context, final Iterator<VariantContext> variants) {
		 Iterator<VariantContextDirectedBreakpoint> it = Iterators.transform(variants, new Function<VariantContext, VariantContextDirectedBreakpoint>() {
			public VariantContextDirectedBreakpoint apply(VariantContext v) {
				return new VariantContextDirectedBreakpointBuilder(context, v).make();
			}
		});
		it = Iterators.filter(it, new Predicate<VariantContextDirectedBreakpoint>() {
			public boolean apply(VariantContextDirectedBreakpoint v) {
				return v.isValid();
			}
		});
		return it;
	}
}
