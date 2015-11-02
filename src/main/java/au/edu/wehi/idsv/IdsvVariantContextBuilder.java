package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.Set;

import com.google.common.collect.Sets;

import au.edu.wehi.idsv.vcf.SvType;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfConstants;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

/**
 * Builder for generating VCF structural variation calls with appropriate attributes
 * @author Daniel Cameron
 *
 */
public class IdsvVariantContextBuilder extends VariantContextBuilder {
	public static final String SOURCE_NAME = "idsv";
	protected final ProcessingContext processContext;
	private final Set<EvidenceSource> sourceSet = Sets.newHashSet();
	public IdsvVariantContextBuilder(ProcessingContext processContext) {
		super();
		this.processContext = processContext;
	}
	public IdsvVariantContextBuilder(ProcessingContext processContext, VariantContext parent) {
		super(parent);
		this.processContext = processContext;
		if (parent instanceof IdsvVariantContext) {
			sourceSet.add(((VariantContextDirectedEvidence)parent).getEvidenceSource());
		}
	}
	public IdsvVariantContextBuilder source(EvidenceSource source) {
		sourceSet.add(source);
		return this;
	}
	private String getBreakendString() {
		return processContext.getConfig().getVariantCalling().placeholderBreakend ? VcfConstants.VCF41BREAKEND_REPLACEMENT : VcfConstants.VCF42BREAKEND;
	}
	public IdsvVariantContextBuilder referenceReads(int[] count) {
		attribute(VcfAttributes.REFERENCE_READ_COUNT.attribute(), count);
		return this;
	}
	public IdsvVariantContextBuilder referenceSpanningPairs(int[] count) {
		attribute(VcfAttributes.REFERENCE_READPAIR_COUNT.attribute(), count);
		return this;
	}
	/**
	 * Sets the variant to the given breakend.
	 * Supporting evidence is not set
	 * @param loc location of breakend
	 * @param untemplatedSequence untemplated breakend sequence 
	 * @return builder
	 */
	public IdsvVariantContextBuilder breakend(BreakendSummary loc, String untemplatedSequence) {
		return dobreak(loc, untemplatedSequence, null);
	}
	/**
	 * Sets the variant to the given breakend.
	 * Supporting evidence is not set
	 * @param loc location of breakend
	 * @param untemplatedSequence untemplated breakend sequence
	 * @param  untemplatedBaseQual base qualities of untemplated breakend sequence 
	 * @return builder
	 */
	public IdsvVariantContextBuilder breakend(BreakendSummary loc, byte[] untemplatedSequence, byte[] untemplatedBaseQual) {
		return dobreak(loc, untemplatedSequence == null ? "" : new String(untemplatedSequence, StandardCharsets.US_ASCII), untemplatedBaseQual);
	}
	/**
	 * Sets the variant to the given breakpoint
	 * Supporting evidence is not set
	 * @param loc location of breakpoint
	 * @param untemplatedSequence untemplated breakpoint sequence
	 * @return builder
	 */
	public IdsvVariantContextBuilder breakpoint(BreakpointSummary loc, String untemplatedSequence) {
		return dobreak(loc, untemplatedSequence, null);
	}
	/**
	 * Sets the variant to the given breakend or breakpoint
	 * @param loc location of break
	 * @param untemplatedSequence untemplated break sequence
	 * @return builder
	 */
	private IdsvVariantContextBuilder dobreak(BreakendSummary loc, String untemplatedSequence, byte[] untemplatedBaseQual) {
		if (loc == null) throw new IllegalArgumentException("loc must be specified");
		if (untemplatedSequence == null) throw new IllegalArgumentException("untemplatedSequence must be specified");
		if (untemplatedBaseQual != null && untemplatedSequence.length() != untemplatedBaseQual.length) {
			throw new IllegalArgumentException("untemplatedBaseQual length must match untemplatedSequence length");
		}
		if (!loc.isValid(processContext.getDictionary())) throw new IllegalArgumentException(String.format("%s is not valid", loc));
		
		String chr = processContext.getDictionary().getSequence(loc.referenceIndex).getSequenceName();
		String ref, alt;
		// call the middle position of inexact intervals
		BreakendSummary bpCall = loc.getCallPosition();
		int callPos = bpCall.start;
		int remoteCallPos = bpCall instanceof BreakpointSummary ? ((BreakpointSummary)bpCall).start2 : 0;
		byte refBase = processContext.getReference().getSubsequenceAt(chr, callPos, callPos).getBases()[0];
		if (!SequenceUtil.isValidBase(refBase)) {
			// Convert ambiguous bases to Ns as VCF does not support full IUPAC notation
			refBase = 'N';
		}
		ref = new String(new byte[] { refBase }, StandardCharsets.US_ASCII);
		if (bpCall instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)loc;
			char remoteBracket = bp.direction2 == BreakendDirection.Forward ? ']' : '[';
			String target = String.format("%s:%d", processContext.getDictionary().getSequence(bp.referenceIndex2).getSequenceName(), remoteCallPos);
			if (bpCall.direction == BreakendDirection.Forward) {
				alt = String.format("%s%s%c%s%c", ref, untemplatedSequence, remoteBracket, target, remoteBracket);
			} else {
				alt = String.format("%c%s%c%s%s", remoteBracket, target, remoteBracket, untemplatedSequence, ref);
			}
		} else {
			if (loc.direction == BreakendDirection.Forward) {
				alt = ref + untemplatedSequence + getBreakendString();
			} else {
				alt = getBreakendString() + untemplatedSequence + ref;
			}
		}
		// populate
		loc(chr, callPos, callPos);
		alleles(ref, alt);
		attribute(VcfSvConstants.SV_TYPE_KEY, SvType.BND.name());
		if (loc.end != loc.start) {
			// Set confidence interval on the call if we don't have an exact breakpoint position
			attribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY, new int[] {loc.start - callPos, loc.end - callPos});
		} else {
			rmAttribute(VcfSvConstants.CONFIDENCE_INTERVAL_START_POSITION_KEY);
		}
		if (loc instanceof BreakpointSummary) {
			BreakpointSummary bp = (BreakpointSummary)loc;
			if (bp.end2 != bp.start2) {
				attribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute(), new int[] { bp.start2 - remoteCallPos, bp.end2 - remoteCallPos});
			} else {
				rmAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute());
			}
		} else {
			rmAttribute(VcfAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute());
		}
		return this;
	}
	/**
	 * sets the phred-scaled quality score for variant
	 * @param phred phred-scaled quality score for variant
	 * @return this
	 */
	public IdsvVariantContextBuilder phredScore(double phred) {
		if (phred < 0 || Double.isInfinite(phred) || Double.isNaN(phred)) {
			log10PError(VariantContext.NO_LOG10_PERROR);
		} else {
			log10PError(phred / -10);
		}
		return this;
	}
	public IdsvVariantContextBuilder attribute(final VcfAttributes key, final Object value) {
		return attribute(key.attribute(), value);
	}
	 
	public IdsvVariantContextBuilder attribute(final String key, final Object value) {
		// don't write Idsv properties unnecessarily
		if (VcfAttributes.getAttributefromKey(key) != null && isNullEmpty(value)) {
			rmAttribute(key);
		} else {
			super.attribute(key, value);
		}
		return this;
    }
	/**
	 * Checks if the given attribute contains non-default information.
	 * @param value attribute value to check
	 * @return true if the given value is equivalent the default values for the containing type.  
	 */
	@SuppressWarnings("rawtypes")
	private static boolean isNullEmpty(Object value) {
		if (value == null) return true;
		if (value instanceof Integer && ((Integer)value) == 0) return true;
		if (value instanceof Long && ((Long)value) == 0L) return true;
		if (value instanceof Float && ((Float)value) == 0f) return true;
		if (value instanceof Double && ((Double)value) == 0d) return true;
		if (value instanceof String && ((String)value).length() == 0) return true;
		if (value instanceof Collection && ((Collection)value).size() == 0) return true;
		if (value instanceof Iterable) {
			for (Object o : (Iterable)value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof Object[]) {
			for (Object o : (Object[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof int[]) {
			for (int o : (int[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof long[]) {
			for (long o : (long[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof float[]) {
			for (float o : (float[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof double[]) {
			for (double o : (double[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof short[]) {
			for (short o : (short[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		} else if (value instanceof byte[]) {
			for (byte o : (byte[])value) {
				if (!isNullEmpty(o)) return false;
			}
			return true;
		}
		return false;
	}
	@Override
	public IdsvVariantContext make() {
        VariantContext underlying = super.make();
        if (underlying.getEnd() < underlying.getStart() && underlying.getEnd() == -1) {
        	throw new IllegalStateException(String.format("Sanity check failure: stop not set for %s", underlying)); 
        }
        EvidenceSource firstSource = sourceSet.isEmpty() ? null : sourceSet.iterator().next();
        return IdsvVariantContext.create(processContext, firstSource, underlying);
	}
}

