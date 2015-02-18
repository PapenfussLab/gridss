package au.edu.wehi.idsv.validation;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import au.edu.wehi.idsv.BreakendAnnotator;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.StructuralVariationCallBuilder;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Annotates breakends based on a reference truth file
 * @author cameron.d
 *
 */
public class TruthAnnotator extends AbstractIterator<VariantContextDirectedEvidence> implements BreakendAnnotator {
	private final ProcessingContext processContext;
	private final List<IdsvVariantContext> truth;
	private final Iterator<VariantContextDirectedEvidence> it;
	public TruthAnnotator(ProcessingContext processContext, Iterator<VariantContextDirectedEvidence> it, File truthVcf) {
		this.processContext = processContext;
		this.truth = loadTruthVcf(processContext, truthVcf);
		this.it = it;
	}
	private static List<IdsvVariantContext> loadTruthVcf(ProcessingContext processContext, File truthVcf) {
		List<IdsvVariantContext> truth = Lists.newArrayList();
		VCFFileReader vcfReader = null;
		CloseableIterator<VariantContext> it = null;
		try {
			vcfReader = new VCFFileReader(truthVcf, false);
			it = vcfReader.iterator();
			while (it.hasNext()) {
				truth.add(IdsvVariantContext.create(processContext, null, it.next()));
			}
		} finally {
			CloserUtil.close(it);
			CloserUtil.close(vcfReader);
		}
		return truth;
	}
	public VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence variant) {
		BreakendSummary variantBreakend = variant.getBreakendSummary();
		BreakpointSummary variantBreakpoint = null;
		if (variantBreakend instanceof BreakpointSummary) {
			variantBreakpoint = (BreakpointSummary)variantBreakend;
			variantBreakend = variantBreakpoint.localBreakend();
		}
		HashSet<String> breakpointHits = Sets.newHashSet();
		HashSet<String> breakendHits = Sets.newHashSet();
		for (IdsvVariantContext truthVariant : truth) {
			if ((truthVariant.hasAttribute(VcfSvConstants.SV_TYPE_KEY) && truthVariant.getAttributeAsString(VcfSvConstants.SV_TYPE_KEY, "").equals("INS")) ||
					(truthVariant.hasAttribute(VcfSvConstants.SV_TYPE_KEY) && truthVariant.getAttributeAsString(VcfSvConstants.SV_TYPE_KEY, "").equals("DEL"))) {
				int svLen = truthVariant.getAttributeAsInt(VcfSvConstants.SV_LENGTH_KEY, 0);
				int untemplatedSequence = 0;
				if (variant instanceof VariantContextDirectedEvidence) {
					untemplatedSequence = ((VariantContextDirectedEvidence)variant).getBreakendSequence().length;
				}
				// two breakpoints: forward & backward
				BreakendSummary truthBreakendStart = new BreakendSummary(truthVariant.getReferenceIndex(), BreakendDirection.Forward, truthVariant.getStart(), truthVariant.getStart());
				int endOffset = Math.max(0, -truthVariant.getAttributeAsInt(VcfSvConstants.SV_LENGTH_KEY, 0));
				endOffset = truthVariant.getStart() + endOffset + 1;
				BreakendSummary truthBreakendEnd = new BreakendSummary(truthVariant.getReferenceIndex(), BreakendDirection.Backward, endOffset, endOffset);
				BreakpointSummary truthBreakpoint = new BreakpointSummary(truthBreakendStart, truthBreakendEnd);
				
				if (matches(truthBreakpoint, variantBreakpoint, svLen, untemplatedSequence)) {
					breakpointHits.add(truthVariant.getID());
				} else if (matches(truthBreakpoint, variantBreakend, svLen, untemplatedSequence)) {
					breakendHits.add(truthVariant.getID());
				}
			} else {
				throw new RuntimeException(String.format("Matching of truth variant at %s:%d not yet implemented.", truthVariant.getChr(), truthVariant.getStart()));
			}
		}
		if (!breakpointHits.isEmpty() || !breakendHits.isEmpty()) {
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, variant);
			builder.attribute(VcfAttributes.TRUTH_MATCHES.attribute(), Lists.newArrayList(breakpointHits));
			builder.attribute(VcfAttributes.TRUTH_MISREALIGN.attribute(), Lists.newArrayList(breakendHits));
			return (VariantContextDirectedEvidence) builder.make();
		}
		return variant;
	}
	public static boolean matches(BreakpointSummary truthBreakpoint, BreakendSummary variantBreakend, int svLen, int untemplatedSequence) {
		return variantBreakend.overlaps(addMargin(truthBreakpoint.localBreakend(), svLen, untemplatedSequence)) ||
				variantBreakend.overlaps(addMargin(truthBreakpoint.remoteBreakend(), svLen, untemplatedSequence));
	}
	public static boolean matches(BreakpointSummary truthBreakpoint, BreakpointSummary variantBreakpoint, int svLen, int untemplatedSequence) {
		truthBreakpoint = new BreakpointSummary(addMargin(truthBreakpoint, svLen, untemplatedSequence), addMargin(truthBreakpoint.remoteBreakend(), svLen, untemplatedSequence));
		return variantBreakpoint.overlaps(truthBreakpoint) || variantBreakpoint.overlaps(truthBreakpoint.remoteBreakpoint());
	}
	private static BreakendSummary addMargin(BreakendSummary loc, int svLen, int untemplatedSequence) {
		int margin = Math.max(0, svLen - untemplatedSequence);
		margin = Math.min(16, margin); // limit to 16bp
		return new BreakendSummary(loc.referenceIndex, loc.direction,
				loc.start - (loc.direction == BreakendDirection.Forward ? 0 : margin),
				loc.end + (loc.direction == BreakendDirection.Backward ? 0 : margin));
	}
	@Override
	protected VariantContextDirectedEvidence computeNext() {
		if (it.hasNext()) return annotate(it.next());
		return endOfData();
	}
}
