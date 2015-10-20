package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.Arrays;
import java.util.List;

import au.edu.wehi.idsv.util.ArrayHelper;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfSvConstants;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private static final Log log = Log.getInstance(StructuralVariationCallBuilder.class);
	private final ProcessingContext processContext;
	private final VariantContextDirectedEvidence parent;
	private final boolean deduplicateEvidence;
	private DirectedBreakpoint bestExactBreakpoint = null;
	private EvidenceIDCollection localEvidence = new EvidenceIDCollection();
	private EvidenceIDCollection assemblyEvidence = new EvidenceIDCollection();
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		this(processContext, parent, true);
	}
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent, boolean deduplicateEvidence) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
		this.deduplicateEvidence = deduplicateEvidence;
	}
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (!isSupportingEvidence(evidence)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		String eid = evidence.getEvidenceID();
		if (localEvidence.get().contains(eid) && deduplicateEvidence) {
			log.debug(String.format("Deduplicating %s from %s", eid, parent.getID()));
			return this;
		}
		localEvidence.categorise(evidence);
		if (evidence instanceof DirectedBreakpoint && evidence.isBreakendExact()) {
			DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
			if (ByBestBreakpointDesc.compare(bp, bestExactBreakpoint) < 0) {
				bestExactBreakpoint = (DirectedBreakpoint)evidence;
			}
		}
		localEvidence.categorise(evidence);
		if (evidence instanceof SAMRecordAssemblyEvidence) {
			if (evidence instanceof RemoteEvidence) {
				assemblyEvidence.addRemote(((SAMRecordAssemblyEvidence)evidence).getEvidenceIDCollection());
			} else {
				assemblyEvidence.add(((SAMRecordAssemblyEvidence)evidence).getEvidenceIDCollection());
			}
		}
		return this;
	}
	private BreakendSummary getBreakendSummaryWithMargin(DirectedEvidence evidence) {
		BreakendSummary bs = evidence.getBreakendSummary();
		bs = processContext.getVariantCallingParameters().withMargin(processContext, bs);
		return bs;
	}
	private boolean isSupportingEvidence(DirectedEvidence evidence) {
		BreakendSummary bs = getBreakendSummaryWithMargin(evidence);
		return parent.getBreakendSummary().overlaps(bs);
	}
	private static float sumPositiveValues(float[] array) {
		// Sorting ascending since all values are positive to minimise error
		Arrays.sort(array);
		// and use doubles to give us extra precision
		double total = 0;
		for (double x : array) {
			total += x;
		}
		return (float)total;
	}
	public VariantContextDirectedEvidence make() {
		int categoryCount = processContext.getCategoryCount();
		//extractAssemblySupport(); // don't extract - as extraction inflates evidence beyond original call
		attribute(VcfAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		attribute(VcfAttributes.BREAKEND_QUAL.attribute(),
				sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKEND_ASSEMBLY_QUAL, categoryCount))
				+ sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKEND_SOFTCLIP_QUAL, categoryCount))
				+ sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKEND_UNMAPPEDMATE_QUAL, categoryCount))
				);
		phredScore(
				sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL, categoryCount))
				+ sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL, categoryCount))
				+ sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKPOINT_READPAIR_QUAL, categoryCount))
				+ sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, categoryCount))
				+ sumPositiveValues(localEvidence.getQual(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL_REMOTE, categoryCount))
				);
		for (VcfAttributes attr : new VcfAttributes[] {
				VcfAttributes.BREAKPOINT_SPLITREAD_COUNT,
				VcfAttributes.BREAKPOINT_READPAIR_COUNT,
				VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE,
				VcfAttributes.BREAKEND_SOFTCLIP_COUNT,
				VcfAttributes.BREAKEND_UNMAPPEDMATE_COUNT
			}) {
			attribute(attr, localEvidence.getCount(attr, categoryCount));
		}
		for (VcfAttributes attr : new VcfAttributes[] {
				VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT,
				VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE,
				VcfAttributes.BREAKEND_ASSEMBLY_COUNT,
			}) {
			attribute(attr, localEvidence.getCount(attr, categoryCount)[0]);
		}
		for (VcfAttributes attr : new VcfAttributes[] {
				VcfAttributes.BREAKPOINT_SPLITREAD_QUAL,
				VcfAttributes.BREAKPOINT_READPAIR_QUAL,
				VcfAttributes.BREAKPOINT_SPLITREAD_QUAL_REMOTE,
				VcfAttributes.BREAKEND_SOFTCLIP_QUAL,
				VcfAttributes.BREAKEND_UNMAPPEDMATE_QUAL
			}) {
			attribute(attr, localEvidence.getQual(attr, categoryCount));
		}
		for (VcfAttributes attr : new VcfAttributes[] {
				VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL,
				VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE,
				VcfAttributes.BREAKEND_ASSEMBLY_QUAL,
			}) {
			attribute(attr, localEvidence.getQual(attr, categoryCount)[0]);
		}
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, assemblyEvidence.getCount(VcfAttributes.BREAKPOINT_READPAIR_COUNT, categoryCount));
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_SPLITREAD_COUNT, ArrayHelper.add(
				assemblyEvidence.getCount(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT, categoryCount),
				assemblyEvidence.getCount(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE, categoryCount)));
		EvidenceIDCollection assemblyUniqueEvidence = new EvidenceIDCollection();
		assemblyUniqueEvidence.add(assemblyEvidence);
		for (String e : localEvidence.get()) {
			assemblyUniqueEvidence.remove(e);
		}
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT, assemblyUniqueEvidence.getCount(VcfAttributes.BREAKPOINT_READPAIR_COUNT, categoryCount));
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT, ArrayHelper.add(
				assemblyUniqueEvidence.getCount(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT, categoryCount),
				assemblyUniqueEvidence.getCount(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE, categoryCount)));
		
		String homo = "";
		if (bestExactBreakpoint != null) {
			breakpoint(bestExactBreakpoint.getBreakendSummary(), bestExactBreakpoint.getUntemplatedSequence());
			homo = bestExactBreakpoint.getHomologySequence();
			rmAttribute(VcfSvConstants.IMPRECISE_KEY);
		} else {
			attribute(VcfSvConstants.IMPRECISE_KEY, true);
		}
		if (homo.length() > 0) {
			attribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY, homo);
			attribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY, homo.length());
		} else {
			rmAttribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY);
			rmAttribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY);
		}

		// id(parent.getID()); // can't change from parent ID as the id is already referenced in the MATEID of the other breakend  
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = applyFilters(variant);
		variant = Models.calculateSomatic(variant);
		//assert(sanitycheck(variant));
		return variant;
	}
	public VariantContextDirectedEvidence applyFilters(VariantContextDirectedEvidence variant) {
		List<VcfFilter> filters = processContext.getVariantCallingParameters().calculateBreakendFilters(variant);
		if (variant instanceof VariantContextDirectedBreakpoint) {
			filters.addAll(processContext.getVariantCallingParameters().calculateBreakpointFilters((VariantContextDirectedBreakpoint)variant));
		}
		if (!filters.isEmpty()) {
			VariantContextBuilder builder = new VariantContextBuilder(variant);
			for (VcfFilter f : filters) {
				builder.filter(f.filter());
			}
			variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, variant.source, builder.make());
		}
		return variant;
	}
	private static Ordering<DirectedBreakpoint> ByBestBreakpointDesc = new Ordering<DirectedBreakpoint>() {
		@Override
		public int compare(DirectedBreakpoint left, DirectedBreakpoint right) {
			return ComparisonChain.start()
					.compareTrueFirst(left.isBreakendExact(), right.isBreakendExact())
					.compareTrueFirst(left instanceof AssemblyEvidence, right instanceof AssemblyEvidence)
					.compare(right.getBreakpointQual(), left.getBreakpointQual())
					.result();
		}
	}.nullsLast();
}
