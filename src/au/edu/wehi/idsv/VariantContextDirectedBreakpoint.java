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
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;

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
	private final byte[] breakendBaseQual;
	//private static Log LOG = Log.getInstance(VariantContextDirectedBreakpoint.class);
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, EvidenceSource source, VariantContext context) {
		this(processContext, source, context, null);
	}
	public VariantContextDirectedBreakpoint(ProcessingContext processContext, EvidenceSource source, VariantContext context, byte[] breakendBaseQual) {
		super(processContext, source, context);
		this.breakendBaseQual = breakendBaseQual;
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
		int rciStart = 0, rciEnd = 0;
		if (hasAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY)) {
			List<Integer> ci = getAttributeAsIntList(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY);
			if (ci.size() == 2) {
				ciStart = ci.get(0);
				ciEnd = ci.get(1);
			} else {
				throw new IllegalStateException(String.format("Error parsing attribute %s of %s. Expected 2 integer values, found %d", VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, super.toString(), ci.size()));
			}
		}
		if (hasAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute())) {
			List<Integer> rci = getAttributeAsIntList(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute());
			if (rci.size() == 2) {
				rciStart = rci.get(0);
				rciEnd = rci.get(1);
			} else {
				throw new IllegalStateException(String.format("Error parsing attribute %s of %s. Expected 2 integer values, found %d", VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute(), super.toString(), rci.size()));
			}
		}
		if (remoteDirection != null) {
			location = new BreakpointSummary(getReferenceIndex(), direction, localPosition + ciStart, localPosition + ciEnd,
					processContext.getDictionary().getSequenceIndex(remoteContig), remoteDirection, remotePosition + rciStart, remotePosition + rciEnd,
					getEvidence());
		} else {
			location = new BreakendSummary(getReferenceIndex(), direction, localPosition + ciStart, localPosition + ciEnd,
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
		return breakendBaseQual;
	}
	@Override
	public boolean isValid() {
		return location != null;
	}
	public String getAssemblerProgram() { return getAttributeAsString(VcfAttributes.ASSEMBLY_PROGRAM.attribute(), null); }
	public String getAssemblyConsensus() { return getAttributeAsString(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), ""); }
	public double getAssemblyQuality() { return getAttributeAsDouble(VcfAttributes.ASSEMBLY_QUALITY.attribute(), 0); }
	public int getAssemblyMaximumSoftClipLength() { return getAttributeAsInt(VcfAttributes.ASSEMBLY_MAX_READ_SOFT_CLIP.attribute(), 0); }
	public int getAssemblyLongestSupportingRead() { return getAttributeAsInt(VcfAttributes.ASSEMBLY_READ_LENGTH.attribute(), 0); }
	public int getTotalReferenceReadDepth() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_READ_COUNT.attribute(), 0, 0); }
	public int getNormalReferenceReadDepth() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_READ_COUNT.attribute(), 1, 0); }
	public int getTumourReferenceReadDepth() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_READ_COUNT.attribute(), 2, 0); }
	public int getTotalReferenceSpanningPairCount() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_SPANNING_READ_PAIR_COUNT.attribute(), 0, 0); }
	public int getNormalReferenceSpanningPairCount() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_SPANNING_READ_PAIR_COUNT.attribute(), 1, 0); }
	public int getTumourReferenceSpanningPairCount() { return getAttributeAsIntListOffset(VcfAttributes.REFERENCE_SPANNING_READ_PAIR_COUNT.attribute(), 2, 0); }
	/**
	 * Returns an iterator containing only the breakend variants from the given iterator
	 * @param context processing context
	 * @param variants input variants
	 * @return breakend variants
	 */
	public static Iterator<VariantContextDirectedBreakpoint> breakendIterator(final ProcessingContext context, final EvidenceSource source, final Iterator<VariantContext> variants) {
		 Iterator<VariantContextDirectedBreakpoint> it = Iterators.transform(variants, new Function<VariantContext, VariantContextDirectedBreakpoint>() {
			public VariantContextDirectedBreakpoint apply(VariantContext v) {
				return new VariantContextDirectedBreakpointBuilder(context, source, v).make();
			}
		});
		it = Iterators.filter(it, new Predicate<VariantContextDirectedBreakpoint>() {
			public boolean apply(VariantContextDirectedBreakpoint v) {
				return v.isValid();
			}
		});
		return it;
	}
	public static Ordering<VariantContextDirectedBreakpoint> ByBreakendStartEnd = new Ordering<VariantContextDirectedBreakpoint>() {
		public int compare(VariantContextDirectedBreakpoint o1, VariantContextDirectedBreakpoint o2) {
			return BreakendSummary.ByStartEnd.compare(o1.getBreakendSummary(), o2.getBreakendSummary());
		  }
	};	
}
