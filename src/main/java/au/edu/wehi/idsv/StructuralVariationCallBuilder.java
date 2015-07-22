package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.util.HashSet;
import java.util.List;

import au.edu.wehi.idsv.sam.SamTags;
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
	private int asCount = 0;
	private int[] scCount = null;
	private int[] rpCount = null;
	private int rasCount = 0;
	private int[] rscCount = null;
	private float totalQual = 0;
	private float asQual = 0;
	private float[] scQual = null;
	private float[] rpQual = null;
	private float rasQual = 0;
	private float[] rscQual;
	private int beasCount = 0;
	private int[] bescCount = null;
	private int[] berpCount = null;
	private float betotalQual = 0;
	private float beasQual = 0;
	private float[] bescQual = null;
	private float[] berpQual = null;
	private HashSet<String> added = new HashSet<String>();
	private EvidenceIDCollection localAssemblyEvidence = new EvidenceIDCollection();
	private EvidenceIDCollection remoteAssemblyEvidence = new EvidenceIDCollection();
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		this(processContext, parent, true);
	}
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent, boolean deduplicateEvidence) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
		this.deduplicateEvidence = deduplicateEvidence;
		scCount = new int[processContext.getCategoryCount()];
		rpCount = new int[processContext.getCategoryCount()];
		rscCount = new int[processContext.getCategoryCount()];
		scQual = new float[processContext.getCategoryCount()];
		rpQual = new float[processContext.getCategoryCount()];
		rscQual = new float[processContext.getCategoryCount()];
		bescCount = new int[processContext.getCategoryCount()];
		berpCount = new int[processContext.getCategoryCount()];
		bescQual = new float[processContext.getCategoryCount()];
		berpQual = new float[processContext.getCategoryCount()];
		added = new HashSet<String>();
	}
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (!isSupportingEvidence(evidence)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		String eid = evidence.getEvidenceID();
		if (added.contains(eid) && deduplicateEvidence) {
			log.debug(String.format("Deduplicating %s from %s", eid, parent.getID()));
			return this;
		}
		added.add(eid);
		if (evidence instanceof DirectedBreakpoint && evidence.isBreakendExact()) {
			DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
			if (ByBestBreakpointDesc.compare(bp, bestExactBreakpoint) < 0) {
				bestExactBreakpoint = (DirectedBreakpoint)evidence;
			}
		}
		if (evidence instanceof SAMRecordAssemblyEvidence) {
			if (evidence instanceof RemoteEvidence) {
				remoteAssemblyEvidence.add(((SAMRecordAssemblyEvidence)evidence).getEvidenceIDCollection());
			} else {
				localAssemblyEvidence.add(((SAMRecordAssemblyEvidence)evidence).getEvidenceIDCollection());
			}
		}
		accumulateBreakendAttributes(evidence);
		if (evidence instanceof DirectedBreakpoint) {
			accumulateBreakpointAttributes((DirectedBreakpoint)evidence);
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
	public VariantContextDirectedEvidence make() {
		//extractAssemblySupport(); // don't extract - as extraction inflates evidence beyond original call
		attribute(VcfAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		phredScore(totalQual);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT.attribute(), asCount); 
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute(), scCount);
		attribute(VcfAttributes.BREAKPOINT_READPAIR_COUNT.attribute(), rpCount);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE.attribute(), rasCount); 
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE.attribute(), rscCount);		
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL.attribute(), asQual); 
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL.attribute(), scQual);
		attribute(VcfAttributes.BREAKPOINT_READPAIR_QUAL.attribute(), rpQual);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE.attribute(), rasQual); 
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL_REMOTE.attribute(), rscQual);
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_COUNT.attribute(), beasCount); 
		attribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT.attribute(), bescCount);
		attribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_COUNT.attribute(), berpCount);
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_QUAL.attribute(), beasQual); 
		attribute(VcfAttributes.BREAKEND_SOFTCLIP_QUAL.attribute(), bescQual);
		attribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_QUAL.attribute(), berpQual);
		attribute(VcfAttributes.BREAKEND_QUAL.attribute(), betotalQual);

		int[] asrp = new int[processContext.getCategoryCount()];
		int[] assr = new int[processContext.getCategoryCount()];
		int[] asrsr = new int[processContext.getCategoryCount()];
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			asrp[i] = localAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_DISCORDANT_PAIR);
			assr[i] = localAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ);
			asrsr[i] = localAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE);
		}
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, asrp);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_SPLITREAD_COUNT, assr);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_REMOTE_SPLITREAD_COUNT, asrsr);
		int[] rasrp = new int[processContext.getCategoryCount()];
		int[] rassr = new int[processContext.getCategoryCount()];
		int[] rasrsr = new int[processContext.getCategoryCount()];
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			rasrp[i] = remoteAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_DISCORDANT_PAIR);
			rassr[i] = remoteAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ);
			rasrsr[i] = remoteAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE);
		}
		attribute(VcfAttributes.BREAKPOINT_REMOTE_ASSEMBLY_READPAIR_COUNT, rasrp);
		attribute(VcfAttributes.BREAKPOINT_REMOTE_ASSEMBLY_SPLITREAD_COUNT, rassr);
		attribute(VcfAttributes.BREAKPOINT_REMOTE_ASSEMBLY_REMOTE_SPLITREAD_COUNT, rasrsr);
		localAssemblyEvidence.removeAll(added);
		// FIXME: remote assembly evidenceIDs are for the other size!
		// Need to switch evidence /1 and /2 for dp and sr with rsr
		remoteAssemblyEvidence.removeAll(added);
		asrp = new int[processContext.getCategoryCount()];
		assr = new int[processContext.getCategoryCount()];
		asrsr = new int[processContext.getCategoryCount()];
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			asrp[i] = localAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_DISCORDANT_PAIR);
			assr[i] = localAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ);
			asrsr[i] = localAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE);
		}
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_READPAIR_CONSCRIPTED_COUNT, asrp);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_SPLITREAD_CONSCRIPTED_COUNT, assr);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_REMOTE_SPLITREAD_CONSCRIPTED_COUNT, asrsr);
		rasrp = new int[processContext.getCategoryCount()];
		rassr = new int[processContext.getCategoryCount()];
		rasrsr = new int[processContext.getCategoryCount()];
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			rasrp[i] = remoteAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_DISCORDANT_PAIR);
			rassr[i] = remoteAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ);
			rasrsr[i] = remoteAssemblyEvidence.getCount(i, SamTags.ASSEMBLY_EVIDENCEID_PREFIX_SPLIT_READ_REMOTE);
		}
		attribute(VcfAttributes.BREAKPOINT_REMOTE_ASSEMBLY_READPAIR_CONSCRIPTED_COUNT, rasrp);
		attribute(VcfAttributes.BREAKPOINT_REMOTE_ASSEMBLY_SPLITREAD_CONSCRIPTED_COUNT, rassr);
		attribute(VcfAttributes.BREAKPOINT_REMOTE_ASSEMBLY_REMOTE_SPLITREAD_CONSCRIPTED_COUNT, rasrsr);
		
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
		variant = calcSpv(variant);
		variant = applyFilters(variant);
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
	private int offsetTN(DirectedEvidence e) {
		EvidenceSource source = e.getEvidenceSource();
		if (source instanceof SAMEvidenceSource) return ((SAMEvidenceSource)source).getSourceCategory();
		return 0;
	}
	private void accumulateBreakendAttributes(DirectedEvidence e) {
		if (e instanceof DirectedBreakpoint) return;
		float qual = e.getBreakendQual();
		int offset = offsetTN(e);
		if (e instanceof AssemblyEvidence) {
			if (e instanceof RemoteEvidence) {
				throw new RuntimeException("Sanity check failure: remote support for only breakend should not be possible");
				//rasCount++;
				//rasQual += qual;
			} else {
				beasCount++;
				beasQual += qual;
			}
		} else if (e instanceof SoftClipEvidence) {
			if (e instanceof RemoteEvidence) {
				throw new RuntimeException("Sanity check failure: remote support for only breakend should not be possible");
				//rscCount[offset]++;
				//rscQual[offset] += qual;
			} else {
				bescCount[offset]++;
				bescQual[offset] += qual;
			}
		} else if (e instanceof NonReferenceReadPair) {
			berpCount[offset]++;
			berpQual[offset] += qual;
		} else {
			throw new RuntimeException("Unknown evidence type " + e.getClass().toString());
		}
		betotalQual += qual;
	}
	private void accumulateBreakpointAttributes(DirectedBreakpoint e) {
		float qual = e.getBreakpointQual();
		int offset = offsetTN(e);
		if (e instanceof AssemblyEvidence) {
			if (e instanceof RemoteEvidence) {
				rasCount++;
				rasQual += qual;
			} else {
				asCount++;
				asQual += qual;
			}
		} else if (e instanceof SoftClipEvidence) {
			if (e instanceof RemoteEvidence) {
				rscCount[offset]++;
				rscQual[offset] += qual;
			} else {
				scCount[offset]++;
				scQual[offset] += qual;
			}
		} else if (e instanceof NonReferenceReadPair) {
			rpCount[offset]++;
			rpQual[offset] += qual;
		} else {
			throw new RuntimeException("Unknown evidence type " + e.getClass().toString());
		}
		totalQual += qual;
	}
	private VariantContextDirectedEvidence calcSpv(VariantContextDirectedEvidence variant) {
		// TODO: somatic p-value should use evidence from both sides of the breakpoint
		float pvalue = Models.somaticPvalue(processContext, variant);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, variant);
		builder.attribute(VcfAttributes.SOMATIC_P_VALUE, pvalue);
		if (pvalue < processContext.getVariantCallingParameters().somaticPvalueThreshold) {
			builder.attribute(VCFConstants.SOMATIC_KEY, true);
		} else {
			builder.rmAttribute(VCFConstants.SOMATIC_KEY);
		}
		return (VariantContextDirectedEvidence)builder.make();
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
