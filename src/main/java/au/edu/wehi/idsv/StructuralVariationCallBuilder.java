package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private static final Log log = Log.getInstance(StructuralVariationCallBuilder.class);
	private final ProcessingContext processContext;
	private final VariantContextDirectedEvidence parent;
	private final Set<String> encounteredEvidenceIDs;
	private final List<DirectedBreakpoint> supportingBreakpoint = new ArrayList<>();
	private final List<DirectedEvidence> supportingBreakend = new ArrayList<>();
	// breakpoint support
	private final List<List<SplitReadEvidence>> supportingSR = new ArrayList<>();
	private final List<List<IndelEvidence>> supportingIndel = new ArrayList<>();
	private final List<List<DiscordantReadPair>> supportingDP = new ArrayList<>();
	private final List<SingleReadEvidence> supportingAS = new ArrayList<>();
	private final List<SingleReadEvidence> supportingRAS = new ArrayList<>();
	private final List<SingleReadEvidence> supportingCAS = new ArrayList<>();
	// breakend support
	private final List<List<SoftClipEvidence>> supportingSC = new ArrayList<>();
	private final List<List<UnmappedMateReadPair>> supportingOEA = new ArrayList<>();
	private final List<SoftClipEvidence> supportingBAS = new ArrayList<>();
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		this(processContext, parent, true);
	}
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent, boolean deduplicateEvidence) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
		this.encounteredEvidenceIDs = deduplicateEvidence ? new HashSet<String>() : null;
		ensureGenotypeBuilders(processContext);
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			supportingSR.add(new ArrayList<>());
			supportingIndel.add(new ArrayList<>());
			supportingDP.add(new ArrayList<>());
			supportingSC.add(new ArrayList<>());
			supportingOEA.add(new ArrayList<>());
		}
	}
	private static int deduplicationMessageCount = 0;
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (!isSupportingEvidence(evidence)) {
			if (isBreakend() && evidence instanceof DirectedBreakpoint) {
				throw new IllegalArgumentException(String.format("Sanity check failure: Breakpoint evidence %s %s should not be assigned to support breakend call at %s",
						evidence.getEvidenceID(),
						evidence.getBreakendSummary(),
						parent.getBreakendSummary()));
			} else {
				throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s %s does not provide support for call at %s",
						evidence.getEvidenceID(),
						evidence.getBreakendSummary(),
						parent.getBreakendSummary()));
			}
		}
		String eid = evidence.getEvidenceID();
		if (encounteredEvidenceIDs != null) {
			if (encounteredEvidenceIDs.contains(eid)) {
				if (deduplicationMessageCount < gridss.Defaults.SUPPRESS_DATA_ERROR_MESSAGES_AFTER) { 
					log.debug(String.format("Deduplicating %s from %s", eid, parent.getID()));
					deduplicationMessageCount++;
					if (deduplicationMessageCount == gridss.Defaults.SUPPRESS_DATA_ERROR_MESSAGES_AFTER) {
						log.debug(String.format("Supressing further deduplication log messages."));
					}
				}
				return this;
			}
			encounteredEvidenceIDs.add(eid);
		}
		if (evidence instanceof DirectedBreakpoint) {
			supportingBreakpoint.add((DirectedBreakpoint)evidence);
		} else {
			supportingBreakend.add(evidence);
		}
		int category = ((SAMEvidenceSource)evidence.getEvidenceSource()).getSourceCategory();
		assert(category < processContext.getCategoryCount());
		if (evidence instanceof DiscordantReadPair) {
			supportingDP.get(category).add((DiscordantReadPair)evidence);
		} else if (evidence instanceof UnmappedMateReadPair) {
			supportingOEA.get(category).add((UnmappedMateReadPair)evidence);
		} else if (evidence instanceof SingleReadEvidence) {
			SingleReadEvidence sre = (SingleReadEvidence) evidence; 
			if (AssemblyAttributes.isAssembly(sre)) {
				AssemblyAttributes attr = new AssemblyAttributes(sre);
				if (evidence instanceof SoftClipEvidence) {
					supportingBAS.add((SoftClipEvidence)sre);
				} else if (!sre.involvesPrimaryReadAlignment()) {
					supportingCAS.add((SingleReadEvidence)sre);
				} else if (sre.getSAMRecord().getSupplementaryAlignmentFlag() || attr.getAssemblyDirection() != sre.getBreakendSummary().direction) {
					supportingRAS.add((SingleReadEvidence)sre);
				} else {
					supportingAS.add((SingleReadEvidence)sre);
				}
			} else if (sre instanceof SoftClipEvidence) {
				supportingSC.get(category).add((SoftClipEvidence)sre);
			} else if (sre instanceof SplitReadEvidence) {
				supportingSR.get(category).add((SplitReadEvidence)sre);
			} else if (sre instanceof IndelEvidence) {
				supportingIndel.get(category).add((IndelEvidence) evidence);
			} else {
				throw new IllegalArgumentException("Unknown evidence type " + evidence.getClass().getName());
			}
		} else {
			throw new IllegalArgumentException("Unknown evidence type " + evidence.getClass().getName());
		}
		return this;
	}
	private static void processAnchor(RangeSet<Integer> anchoredBases, SAMRecord record) {
		if (record.getReadUnmappedFlag()) return;
		if (CigarUtil.widthOfImprecision(record.getCigar()) > 0) return; // unanchored
		int pos = record.getAlignmentStart();
		for (CigarElement ce : record.getCigar().getCigarElements()) {
			switch (ce.getOperator()) {
				case M:
				case EQ:
				case X:
					anchoredBases.add(Range.closedOpen(pos, pos + ce.getLength()));
					break;
				default:
					break;
			}
			if (ce.getOperator().consumesReferenceBases()) {
				pos += ce.getLength();
			}
		}
	}
	/**
	 * Creates a cigar for alignment to given breakend anchored bases, with an XNX alignment to the breakend interval
	 * @param anchors anchored bases
	 * @param bp breakend
	 * @return Cigar representing breakend with anchored bases
	 */
	private static Cigar makeCigar(RangeSet<Integer> anchors, BreakendSummary bp) {
		if (bp.direction == BreakendDirection.Forward) {
			anchors = anchors.subRangeSet(Range.closedOpen(Integer.MIN_VALUE, bp.end + 1));
		} else {
			anchors = anchors.subRangeSet(Range.closedOpen(bp.start, Integer.MAX_VALUE));
		}
		anchors = TreeRangeSet.create(anchors);
		anchors.add(Range.closed(bp.start, bp.end + 1));
		List<CigarElement> cigar = new ArrayList<CigarElement>();
		int lastEndPos = anchors.span().lowerEndpoint();
		for (Range<Integer> span : anchors.asRanges()) {
			cigar.add(new CigarElement(span.lowerEndpoint() - lastEndPos, CigarOperator.D));
			cigar.add(new CigarElement(span.upperEndpoint() - span.lowerEndpoint(), CigarOperator.M));
			lastEndPos = span.upperEndpoint();
		}
		cigar.remove(0); //cigar = CigarUtil.clean(cigar);
		// replace the placeholder homology alignment with XNX notation
		List<CigarElement> breakcigar = bp.getCigarRepresentation();
		int width = bp.end - bp.start + 1;
		if (bp.direction == BreakendDirection.Forward) {
			// trim 
			CigarElement e = cigar.remove(cigar.size() - 1);
			if (e.getLength() > width) {
				cigar.add(new CigarElement(e.getLength() - width, e.getOperator()));
			}
			cigar.addAll(breakcigar);
		} else {
			CigarElement e = cigar.remove(0);
			if (e.getLength() > width) {
				cigar.add(0, new CigarElement(e.getLength() - width, e.getOperator()));
			}
			cigar.addAll(0, breakcigar);
		}
		return new Cigar(cigar);
	}
	private BreakendSummary getBreakendSummaryWithMargin(DirectedEvidence evidence) {
		BreakendSummary bs = evidence.getBreakendSummary();
		bs = processContext.getVariantCallingParameters().withMargin(bs);
		return bs;
	}
	private boolean isSupportingEvidence(DirectedEvidence evidence) {
		BreakendSummary bs = getBreakendSummaryWithMargin(evidence);
		return parent.getBreakendSummary().overlaps(bs) &&
				// breakpoint evidence should not be assigned to breakend calls
				!(isBreakend() && evidence instanceof DirectedBreakpoint);
	}
	private boolean isBreakend() {
		return !(parent.getBreakendSummary() instanceof BreakpointSummary);
	}
	public VariantContextDirectedEvidence make() {
		attribute(VcfInfoAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		double beQual = supportingBAS.stream().mapToDouble(e -> e.getBreakendQual()).sum() 
				+ supportingSC.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakendQual()).sum()
				+ supportingOEA.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakendQual()).sum();
		double bpQual = supportingAS.stream().mapToDouble(e -> ((DirectedBreakpoint)e).getBreakpointQual()).sum()
				+ supportingRAS.stream().mapToDouble(e -> ((DirectedBreakpoint)e).getBreakpointQual()).sum()
				+ supportingCAS.stream().mapToDouble(e -> ((DirectedBreakpoint)e).getBreakpointQual()).sum()
				+ supportingSR.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakpointQual()).sum()
				+ supportingIndel.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakpointQual()).sum()
				+ supportingDP.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakpointQual()).sum();
		attribute(VcfInfoAttributes.BREAKEND_QUAL.attribute(), beQual);
		phredScore(isBreakend() ? beQual : bpQual);
		// Count
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT, supportingAS.size());
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE, supportingRAS.size());
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_COMPOUND, supportingCAS.size());
		attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_COUNT, supportingBAS.size());
		sumIntAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT, VcfFormatAttributes.BREAKPOINT_SPLITREAD_COUNT, supportingSR, e -> 1);
		sumIntAttr(VcfInfoAttributes.BREAKPOINT_INDEL_COUNT, VcfFormatAttributes.BREAKPOINT_INDEL_COUNT, supportingIndel, e -> 1);
		sumIntAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_READPAIR_COUNT, supportingDP, e -> 1);
		sumIntAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT, VcfFormatAttributes.BREAKEND_SOFTCLIP_COUNT, supportingSC, e -> 1);
		sumIntAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_COUNT, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_COUNT, supportingOEA, e -> 1);
		// Qual
		double[] asq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_QUAL, supportingAS);
		double[] rasq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, supportingRAS);
		double[] casq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_COMPOUND, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_QUAL_COMPOUND, supportingCAS);
		double[] basq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKEND_ASSEMBLY_QUAL, VcfFormatAttributes.BREAKEND_ASSEMBLY_QUAL, supportingBAS);
		double[] srq = sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_QUAL, VcfFormatAttributes.BREAKPOINT_SPLITREAD_QUAL, supportingSR, e -> e.getBreakpointQual());
		double[] iq = sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_INDEL_QUAL, VcfFormatAttributes.BREAKPOINT_INDEL_QUAL, supportingIndel, e -> e.getBreakpointQual());
		double[] rpq = sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_QUAL, VcfFormatAttributes.BREAKPOINT_READPAIR_QUAL, supportingDP, e -> e.getBreakpointQual());
		double[] scq = sumDoubleAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_QUAL, VcfFormatAttributes.BREAKEND_SOFTCLIP_QUAL, supportingSC, e -> e.getBreakendQual());
		double[] umq = sumDoubleAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_QUAL, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_QUAL, supportingOEA, e -> e.getBreakendQual());
		//sumAttr(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT, fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT);
		//sumAttr(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READ_COUNT, fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT);
		
		// Assembly breakdown
		int[] asr = new int[processContext.getCategoryCount()];
		int[] asrp = new int[processContext.getCategoryCount()];
		int[] basr = new int[processContext.getCategoryCount()];
		int[] basrp = new int[processContext.getCategoryCount()];
		int[] supportingBreakpointFragments = new int[processContext.getCategoryCount()];
		int[] supportingBreakendFragments = new int[processContext.getCategoryCount()];
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			int category = i;
			asr[category] = Stream.concat(Stream.concat(supportingAS.stream(), supportingRAS.stream()), supportingCAS.stream())
					.mapToInt(ass -> ass.getCategorySupportBreakdown().get(category) ? new AssemblyAttributes(ass.getSAMRecord()).getAssemblySupportCountSoftClip(category) : 0)
					.sum();
			asrp[category] = Stream.concat(Stream.concat(supportingAS.stream(), supportingRAS.stream()), supportingCAS.stream())
					.mapToInt(ass -> ass.getCategorySupportBreakdown().get(category) ? new AssemblyAttributes(ass.getSAMRecord()).getAssemblySupportCountReadPair(category) : 0)
					.sum();
			basr[category] = supportingBAS.stream()
					.mapToInt(ass -> ass.getCategorySupportBreakdown().get(category) ? new AssemblyAttributes(ass.getSAMRecord()).getAssemblySupportCountSoftClip(category) : 0)
					.sum();
			basrp[category] = supportingBAS.stream()
					.mapToInt(ass -> ass.getCategorySupportBreakdown().get(category) ? new AssemblyAttributes(ass.getSAMRecord()).getAssemblySupportCountReadPair(category) : 0)
					.sum();
			Set<String> bpfrags = supportingBreakpoint.stream()
					.map(e -> e.getOriginatingFragmentID(category))
					.flatMap(x -> x.stream())
					.collect(Collectors.toSet());
			supportingBreakpointFragments[category] = bpfrags.size();
			Set<String> befrags = supportingBreakend.stream()
					.map(e -> e.getOriginatingFragmentID(category))
					.flatMap(x -> x.stream())
					.collect(Collectors.toSet());
			befrags.removeAll(bpfrags);
			supportingBreakendFragments[category] = befrags.size(); 
					
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT.attribute(), asr[category]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT.attribute(), asrp[category]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_ASSEMBLY_READ_COUNT.attribute(), basr[category]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_ASSEMBLY_READPAIR_COUNT.attribute(), basrp[category]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_QUAL.attribute(), srq[i] + iq[i] + rpq[i] + asq[i] + rasq[i] + casq[i]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_QUAL.attribute(), scq[i] + umq[i] + basq[i]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_VARIANT_FRAGMENTS.attribute(), supportingBreakpointFragments[category]);
			genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_VARIANT_FRAGMENTS.attribute(), supportingBreakendFragments[category]);
		}
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT.attribute(), IntStream.of(asr).sum());
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT.attribute(), IntStream.of(asrp).sum());
		attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_READ_COUNT.attribute(), IntStream.of(basr).sum());
		attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_READPAIR_COUNT.attribute(), IntStream.of(basrp).sum());
		attribute(VcfInfoAttributes.BREAKPOINT_VARIANT_FRAGMENTS.attribute(), IntStream.of(supportingBreakpointFragments).sum());
		attribute(VcfInfoAttributes.BREAKEND_VARIANT_FRAGMENTS.attribute(), IntStream.of(supportingBreakendFragments).sum());
		
		if (supportingBreakpoint.size() > 0) {
			int totalSupport = supportingBreakpoint.stream().mapToInt(e -> e.constituentReads()).sum();
			if (totalSupport != 0) {
				attribute(VcfInfoAttributes.STRAND_BIAS.attribute(), supportingBreakpoint.stream().mapToDouble(e -> e.getStrandBias() * e.constituentReads()).sum() / totalSupport);
			}
		} else {
			int totalSupport = supportingBreakend.stream().mapToInt(e -> e.constituentReads()).sum();
			if (totalSupport != 0) {
				attribute(VcfInfoAttributes.STRAND_BIAS.attribute(), supportingBreakend.stream().mapToDouble(e -> e.getStrandBias() * e.constituentReads()).sum() / totalSupport);
			}
		}
		List<String> breakendIds = Stream.concat(Stream.concat(Stream.concat(supportingAS.stream(), supportingRAS.stream()), supportingCAS.stream()), supportingBAS.stream())
				.map(e -> e.getSAMRecord().getReadName())
				.distinct()
				.sorted() // ensure deterministic output
				.collect(Collectors.toList());
		if (breakendIds.size() > 0) {
			attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID, breakendIds);
		} else {
			rmAttribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID.attribute());
		}
		
		String untemplated = parent.getBreakpointSequenceString();
		String homo = "";
		BreakendSummary nominalPosition = parent.getBreakendSummary();
		
		if (isBreakend()) {
			DirectedEvidence bestBreakend = supportingBreakend.stream()
					.sorted(ByBestBreakendDesc)
					.findFirst().orElse(null);
			if (bestBreakend != null && bestBreakend.isBreakendExact()) {
				untemplated = new String(bestBreakend.getBreakendSequence());
				nominalPosition = bestBreakend.getBreakendSummary();
				breakend(bestBreakend.getBreakendSummary(), untemplated);
				rmAttribute(VcfSvConstants.IMPRECISE_KEY);
			} else {
				attribute(VcfSvConstants.IMPRECISE_KEY, true);
			}
		} else {
			DirectedBreakpoint bestBreakpoint = supportingBreakpoint.stream()
					.sorted(ByBestBreakpointDesc)
					.findFirst().orElse(null);
			if (bestBreakpoint != null && bestBreakpoint.isBreakendExact()) {
				untemplated = bestBreakpoint.getUntemplatedSequence();
				nominalPosition = bestBreakpoint.getBreakendSummary();
				breakpoint(bestBreakpoint.getBreakendSummary(), untemplated);
				homo = bestBreakpoint.getHomologySequence();
				rmAttribute(VcfSvConstants.IMPRECISE_KEY);
			} else {
				attribute(VcfSvConstants.IMPRECISE_KEY, true);
			}
		}
		if (homo.length() > 0) {
			attribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY, homo);
			attribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY, homo.length());
		} else {
			rmAttribute(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY);
			rmAttribute(VcfSvConstants.HOMOLOGY_LENGTH_KEY);
		}
		RangeSet<Integer> allAnchoredBases = TreeRangeSet.create();
		// TODO: per sample support cigar
		supportingSR.stream().flatMap(l -> l.stream()).forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingIndel.stream().flatMap(l -> l.stream()).forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingDP.stream().flatMap(l -> l.stream()).forEach(e -> processAnchor(allAnchoredBases, e.getLocalledMappedRead()));
		supportingAS.stream().forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingRAS.stream().forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingCAS.stream().forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		attribute(VcfInfoAttributes.SUPPORT_CIGAR, makeCigar(allAnchoredBases, nominalPosition).toString());
		
		// id(parent.getID()); // can't change from parent ID as the id is already referenced in the MATEID of the other breakend  
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = applyFilters(variant);
		//variant = Models.calculateSomatic(variant);
		return variant;
	}
	/*
	 * To get per-sample assembly scores, we pro-rata the assembly score by the contribution from each category.
	 * 
	 * Note that the sum of the contributing breakend scores is not necessarily the same as the assembly breakpoint score.
	 * This means that each assembly must be individually pro-rataed. 
	 */
	private <T extends SingleReadEvidence> double[] prorataAssemblyQualBreakdown(VcfInfoAttributes infoAttr, VcfFormatAttributes attr, List<T> assemblies) {
		double totalAssQual = 0;
		double[] prorata = new double[processContext.getCategoryCount()];
		for (SingleReadEvidence ass : assemblies) {
			List<Boolean> supportingCategories = ass.getCategorySupportBreakdown();
			AssemblyAttributes aa = new AssemblyAttributes(ass.getSAMRecord());
			double assQual = (ass instanceof DirectedBreakpoint) ? ((DirectedBreakpoint)ass).getBreakpointQual() : ass.getBreakendQual();
			double beQual = aa.getAssemblySupportSoftClipQualityScore(supportingCategories) + aa.getAssemblySupportReadPairQualityScore(supportingCategories);
			for (int i = 0; i < prorata.length; i++) {
				if (supportingCategories.get(i)) {
					prorata[i] += assQual * (aa.getAssemblySupportSoftClipQualityScore(i) + aa.getAssemblySupportReadPairQualityScore(i)) / beQual;
				}
			}
			totalAssQual += assQual;
		}
		for (int i = 0; i < prorata.length; i++) {
			genotypeBuilder.get(i).attribute(attr.attribute(), prorata[i]);
		}
		attribute(infoAttr.attribute(), totalAssQual);
		return prorata;
	}
	private <T> void sumIntAttr(
			VcfInfoAttributes infoAttr,
			VcfFormatAttributes formatAttr,
			List<List<T>> support,
			ToIntFunction<T> f) {
		assert(support.size() == processContext.getCategoryCount());
		int sum = 0;
		for (int i = 0; i < support.size(); i++) {
			int value = support.get(i).stream().mapToInt(f).sum();
			genotypeBuilder.get(i).attribute(formatAttr.attribute(), value);
			sum += value;
		}
		attribute(infoAttr, sum);
	}
	private <T extends DirectedEvidence> double[] sumDoubleAttr(
			VcfInfoAttributes infoAttr,
			VcfFormatAttributes formatAttr,
			List<List<T>> support,
			ToDoubleFunction<T> f) {
		assert(support.size() == processContext.getCategoryCount());
		double sum = 0;
		double[] result = new double[processContext.getCategoryCount()]; 
		for (int i = 0; i < support.size(); i++) {
			double value = support.get(i).stream().mapToDouble(f).sum();
			genotypeBuilder.get(i).attribute(formatAttr.attribute(), value);
			sum += value;
			result[i] = value;
		}
		attribute(infoAttr, sum);
		return result;
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
	private static Ordering<DirectedEvidence> ByBestBreakendDesc = new Ordering<DirectedEvidence>() {
		@Override
		public int compare(DirectedEvidence left, DirectedEvidence right) {
			return ComparisonChain.start()
					.compareTrueFirst(left.isBreakendExact(), right.isBreakendExact())
					.compareTrueFirst(AssemblyAttributes.isAssembly(left), AssemblyAttributes.isAssembly(right))
					.compare(right.getBreakendQual(), left.getBreakendQual()) // desc
					.compare(
							right.getBreakendSequence() == null ? -1 : right.getBreakendSequence().length,
							left.getBreakendSequence() == null ? -1 : left.getBreakendSequence().length) // desc
					.compare(left.getEvidenceID(), right.getEvidenceID())
					.result();
		}
	}.nullsLast();
	private static Ordering<DirectedBreakpoint> ByBestBreakpointDesc = new Ordering<DirectedBreakpoint>() {
		@Override
		public int compare(DirectedBreakpoint left, DirectedBreakpoint right) {
			return ComparisonChain.start()
					.compareTrueFirst(left.isBreakendExact(), right.isBreakendExact())
					.compareTrueFirst(AssemblyAttributes.isAssembly(left), AssemblyAttributes.isAssembly(right))
					.compare(right.getBreakpointQual(), left.getBreakpointQual()) // desc
					// Take the one that aligned the most bases (more likely to remove REF FP calls)
					.compare(left.getUntemplatedSequence().length(), right.getUntemplatedSequence().length())
					// Ensure both sides of the breakpoint choose the same evidence in case of a tie 
					.compare(left.getEvidenceID(), right.getEvidenceID())
					.result();
		}
	}.nullsLast();
}
