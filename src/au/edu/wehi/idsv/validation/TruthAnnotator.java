package au.edu.wehi.idsv.validation;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import au.edu.wehi.idsv.BreakendAnnotator;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.StructuralVariationCallBuilder;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

/**
 * Annotates breakends based on a reference truth file
 * @author cameron.d
 *
 */
public class TruthAnnotator implements BreakendAnnotator {
	private final ProcessingContext processContext;
	private final List<IdsvVariantContext> truth;
	public TruthAnnotator(ProcessingContext processContext, File truthVcf) {
		this.processContext = processContext;
		this.truth = loadTruthVcf(processContext, truthVcf);
	}
	private static List<IdsvVariantContext> loadTruthVcf(ProcessingContext processContext, File truthVcf) {
		List<IdsvVariantContext> truth = Lists.newArrayList();
		VCFFileReader vcfReader = null;
		CloseableIterator<VariantContext> it = null;
		try {
			vcfReader = new VCFFileReader(truthVcf);
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
	@Override
	public VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence variant) {
		BreakendSummary truthBreakend = variant.getBreakendSummary();
		BreakpointSummary truthBreakpoint = null;
		if (truthBreakend instanceof BreakpointSummary) {
			truthBreakpoint = (BreakpointSummary)truthBreakend;
			truthBreakend = truthBreakpoint.localBreakend();
		}
		HashSet<String> breakpointHits = Sets.newHashSet();
		HashSet<String> breakendHits = Sets.newHashSet();
		for (IdsvVariantContext truthVariant : truth) {
			if ((truthVariant.hasAttribute(VcfSvConstants.SV_TYPE_KEY) && truthVariant.getAttributeAsString(VcfSvConstants.SV_TYPE_KEY, "").equals("INS")) ||
					(truthVariant.hasAttribute(VcfSvConstants.SV_TYPE_KEY) && truthVariant.getAttributeAsString(VcfSvConstants.SV_TYPE_KEY, "").equals("DEL"))) {
				// two breakpoints: forward & backward
				BreakendSummary start = new BreakendSummary(truthVariant.getReferenceIndex(), BreakendDirection.Forward, truthVariant.getStart(), truthVariant.getStart());
				int endOffset = Math.max(0, truthVariant.getAttributeAsInt(VcfSvConstants.SV_LENGTH_KEY, 0));
				endOffset = truthVariant.getStart() + endOffset + 1;
				BreakendSummary end = new BreakendSummary(truthVariant.getReferenceIndex(), BreakendDirection.Backward, endOffset, endOffset);
				BreakpointSummary bp = new BreakpointSummary(start, end);
				if (variant.getBreakendSummary().overlaps(bp) || variant.getBreakendSummary().overlaps(bp.remoteBreakpoint())) {
					breakpointHits.add(truthVariant.getID());
				} else if (variant.getBreakendSummary().overlaps(bp.localBreakend()) || variant.getBreakendSummary().overlaps(bp.remoteBreakend())) {
					breakendHits.add(truthVariant.getID());
				}
			} else {
				throw new RuntimeException(String.format("Matching of truth variant at %s:%d not yet implemented.", truthVariant.getChr(), truthVariant.getStart()));
			}
		}
		if (!breakpointHits.isEmpty() || !breakendHits.isEmpty()) {
			StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(processContext, variant);
			builder.attribute("REF_MATCHES", Lists.newArrayList(breakpointHits));
			builder.attribute("REF_MISREALIGN", Lists.newArrayList(breakendHits));
			return builder.make();
		}
		return variant;
	}
}
