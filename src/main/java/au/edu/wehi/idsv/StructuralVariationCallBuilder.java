package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.vcf.VcfFilter;
import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.collect.*;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private static final Log log = Log.getInstance(StructuralVariationCallBuilder.class);
	private final ProcessingContext processContext;
	private final CalledBreakpointPositionLookup calledBreakpointLookup;
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
	private boolean updateReadInformation = true;
	private boolean updateAssemblyInformation= true;
	public StructuralVariationCallBuilder(ProcessingContext processContext, CalledBreakpointPositionLookup lookup, VariantContextDirectedEvidence parent) {
		this(processContext, lookup, parent, true);
	}
	public StructuralVariationCallBuilder(ProcessingContext processContext, CalledBreakpointPositionLookup calledBreakpointLookup, VariantContextDirectedEvidence parent, boolean deduplicateEvidence) {
		super(processContext, parent);
		this.calledBreakpointLookup = calledBreakpointLookup;
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
		Map<SingleReadEvidence, AssemblyAttributes> aaLookup = new HashMap<>();
		Stream.of(supportingAS.stream(), supportingRAS.stream(), supportingCAS.stream(), supportingBAS.stream())
			.flatMap(x -> x)
			.forEach(ass -> {
				AssemblyAttributes aa =  new AssemblyAttributes(ass.getSAMRecord());
				aaLookup.put(ass, aa);
			});
		if (isUpdateVariantQualityScore()) {
			updateVariantQualityScoreAttributes(aaLookup);
		}
		updateVariantSupportCounts(aaLookup);
		if (processContext.getConfig().getVariantCalling().includeSupportingReadNames) {
			updateSupportingReadNames();
		}
		updateStrandBias(aaLookup);
		updateAssemblySupportTracking();
		updateNominalCallPosition();
		
		// id(parent.getID()); // can't change from parent ID as the id is already referenced in the MATEID of the other breakend  
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext.getDictionary(), null, super.make());
		variant = applyFilters(variant);
		//variant = Models.calculateSomatic(variant);
		return variant;
	}

	private void updateNominalCallPosition() {
		if (!isUpdateAssemblyInformation()) {
			// nominal positions are usually from assembly so if we're not updating assemblies then we don't change the call
			return;
		}
		BreakendSummary nominalPosition = parent.getBreakendSummary();
		String untemplated = "";// = parent.getBreakpointSequenceString();
		String homo = parent.getAttributeAsString(VcfSvConstants.HOMOLOGY_SEQUENCE_KEY, "");
		if (isBreakend()) {
			DirectedEvidence bestBreakend = supportingBreakend.stream()
					.sorted(ByBestBreakendDesc)
					.findFirst().orElse(null);
			if (bestBreakend != null && bestBreakend.isBreakendExact()) {
				untemplated = new String(bestBreakend.getBreakendSequence());
				nominalPosition = bestBreakend.getBreakendSummary();
				breakend(nominalPosition, untemplated);
				rmAttribute(VcfSvConstants.IMPRECISE_KEY);
			} else {
				if (isUpdateAssemblyInformation() && isUpdateReadInformation()) {
					attribute(VcfSvConstants.IMPRECISE_KEY, true);
				}
			}
		} else {
			String event = parent.getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, null);
			CalledBreakpointPositionLookup.NominalPosition np = calledBreakpointLookup.removeUpper(event);
			boolean isExact;
			if (np == null) {
				DirectedBreakpoint bestBreakpoint = supportingBreakpoint.stream()
						.sorted(ByBestBreakpointDesc)
						.findFirst().orElse(null);
				if (bestBreakpoint != null && bestBreakpoint.isBreakendExact()) {
					untemplated = bestBreakpoint.getUntemplatedSequence();
					nominalPosition = bestBreakpoint.getBreakendSummary().centreAligned();
					breakpoint((BreakpointSummary) nominalPosition, untemplated);
					homo = bestBreakpoint.getHomologySequence();
					isExact = true;
					rmAttribute(VcfSvConstants.IMPRECISE_KEY);
				} else {
					isExact = false;
					if (isUpdateAssemblyInformation() && isUpdateReadInformation()) {
						attribute(VcfSvConstants.IMPRECISE_KEY, true);
					}
				}
				np = new CalledBreakpointPositionLookup.NominalPosition((BreakpointSummary)nominalPosition, untemplated, homo, isExact);
				calledBreakpointLookup.addLower(event, np);
				// check against parent position as we could have moved
				if (parent.getBreakendSummary() instanceof BreakpointSummary && ((BreakpointSummary)parent.getBreakendSummary()).isHighBreakend()) {
					log.info("CalledBreakpointLookup entry missing for " + event);
				}
			}
			breakpoint(np.nominalPosition, np.insertedSequenced);
			homo = np.homologySequence;
			if (np.isExact) {
				rmAttribute(VcfSvConstants.IMPRECISE_KEY);
			} else {
				if (isUpdateAssemblyInformation() && isUpdateReadInformation()) {
					attribute(VcfSvConstants.IMPRECISE_KEY, true);
				}
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
		supportingAS.stream().filter(e -> !shouldfilterAssemblyFromSupportInterval(e)).forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingRAS.stream().filter(e -> !shouldfilterAssemblyFromSupportInterval(e)).forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingCAS.stream().filter(e -> !shouldfilterAssemblyFromSupportInterval(e)).forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		attribute(VcfInfoAttributes.SUPPORT_CIGAR, makeCigar(allAnchoredBases, nominalPosition).toString());
		RangeSet<Integer> assemblyAnchoredBases = TreeRangeSet.create();
		supportingAS.stream().filter(e -> !shouldfilterAssemblyFromSupportInterval(e)).forEach(e -> processAnchor(allAnchoredBases, e.getSAMRecord()));
		supportingRAS.stream().filter(e -> !shouldfilterAssemblyFromSupportInterval(e)).forEach(e -> processAnchor(assemblyAnchoredBases, e.getSAMRecord()));
		supportingCAS.stream().filter(e -> !shouldfilterAssemblyFromSupportInterval(e)).forEach(e -> processAnchor(assemblyAnchoredBases, e.getSAMRecord()));
		attribute(VcfInfoAttributes.ASSEMBLY_SUPPORT_CIGAR, makeCigar(assemblyAnchoredBases, nominalPosition).toString());
	}

	/**
	 * https://github.com/PapenfussLab/gridss/issues/213
	 * assemblies can assemble fragment that do not actually support the breakpoint
	 * (e.g. nearby reads noise soft-clip). This means we can't use local assemblies
	 * in anchor support interval calculation.
	 * Note: this includes both sides of indel-spanning assemblies as well.
	 */
	private static boolean shouldfilterAssemblyFromSupportInterval(SingleReadEvidence e) {
		// TODO: are both sides of indel assemblies
		return !e.getSAMRecord().isSecondaryOrSupplementary();
	}

	private void updateAssemblySupportTracking() {
		if (isUpdateAssemblyInformation()) {
			List<SingleReadEvidence> suportingAssemblies = new ArrayList<>();
			suportingAssemblies.addAll(supportingAS);
			suportingAssemblies.addAll(supportingRAS);
			suportingAssemblies.addAll(supportingCAS);
			suportingAssemblies.addAll(supportingBAS);
			suportingAssemblies.sort(Comparator.comparing(o -> o.getSAMRecord().getReadName()));
			// Track long read support as well
			for (List<SplitReadEvidence> list : supportingSR) {
				for (SingleReadEvidence read : list) {
					if (read.getEvidenceSource().isLongReadLibrary()) {
						suportingAssemblies.add(read);
					}
				}
			}
			if (suportingAssemblies.size() > 0) {
				attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID, suportingAssemblies.stream()
						.map(o -> o.getSAMRecord().getReadName())
						.collect(Collectors.toList()));
				attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID_LOCAL_CONTIG_OFFSET, suportingAssemblies.stream()
						.map(o -> o.getLocalChimericAlignmentReadOffset())
						.collect(Collectors.toList()));
				attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID_REMOTE_CONTIG_OFFSET, suportingAssemblies.stream()
						.map(o -> (o instanceof SplitReadEvidence) ? ((SplitReadEvidence) o).getRemoteChimericAlignmentReadOffset() :
								(o instanceof IndelEvidence ? o.getLocalChimericAlignmentReadOffset() : -1))
						.collect(Collectors.toList()));
			} else {
				rmAttribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID.attribute());
				rmAttribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID_LOCAL_CONTIG_OFFSET.attribute());
				rmAttribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID_REMOTE_CONTIG_OFFSET.attribute());
			}
		}
	}

	private void updateSupportingReadNames() {
		if (isUpdateAssemblyInformation() & isUpdateReadInformation()) {
			List<String> supportingBreakpointNames = Streams.concat(
					supportingIndel.stream().flatMap(l -> l.stream()).map(e -> e.getSAMRecord().getReadName()),
					supportingSR.stream().flatMap(l -> l.stream()).map(e -> e.getSAMRecord().getReadName()),
					supportingDP.stream().flatMap(l -> l.stream()).map(e -> e.getLocalledMappedRead().getReadName()))
					.distinct()
					.sorted()
					.collect(Collectors.toList());
			List<String> supportingBreakendNames = Streams.concat(
					supportingSC.stream().flatMap(l -> l.stream()).map(e -> e.getSAMRecord().getReadName()),
					supportingOEA.stream().flatMap(l -> l.stream()).map(e -> e.getLocalledMappedRead().getReadName()))
					.distinct()
					.sorted()
					.collect(Collectors.toList());
			attribute(VcfInfoAttributes.SUPPORTING_BREAKPOINT_READ_NAMES.attribute(), supportingBreakpointNames);
			attribute(VcfInfoAttributes.SUPPORTING_BREAKEND_READ_NAMES.attribute(), supportingBreakendNames);
		}
	}

	private void updateStrandBias(Map<SingleReadEvidence, AssemblyAttributes> aaLookup) {
		if (isUpdateReadInformation()) {
			// Calculate strand bias purely from direct read support
			int reads = 0;
			float strandReads = 0;
			for (DirectedEvidence de : Iterables.concat(supportingBreakpoint, supportingBreakend)) {
				if (de instanceof SingleReadEvidence) {
					SingleReadEvidence e = (SingleReadEvidence) de;
					if (AssemblyAttributes.isAssembly(e)) {
						AssemblyAttributes aa = aaLookup.get(de);
						int asmReads = aa.getSupportingReadCount(null, null, ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read), null, processContext);
						reads += asmReads;
						strandReads += aa.getStrandBias() * asmReads;
					} else {
						reads++;
						strandReads += e.getStrandBias();
					}
				}
			}
			if (reads > 0) {
				attribute(VcfInfoAttributes.STRAND_BIAS.attribute(), strandReads / reads);
			}
		}
	}

	private void updateVariantSupportCounts(Map<SingleReadEvidence, AssemblyAttributes> aaLookup) {
		if (isUpdateAssemblyInformation() & isUpdateReadInformation()) {
			// Non-supporting assemblies
			Set<String> assNames = supportingAS.stream()
					.map(e -> e.getAssociatedAssemblyName())
					.collect(Collectors.toSet());
			sumIntAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_READPAIR_COUNT,
					supportingDP, e -> (e.getAssociatedAssemblyName() != null && !assNames.contains(e.getAssociatedAssemblyName())) ? 1 : 0);
			sumIntAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_SPLITREAD_COUNT, VcfFormatAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_SPLITREAD_COUNT,
					supportingSR, e -> (e.getAssociatedAssemblyName() != null && !assNames.contains(e.getAssociatedAssemblyName())) ? 1 : 0);
			sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_READPAIR_QUAL, VcfFormatAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_READPAIR_QUAL,
					supportingDP, e -> (e.getAssociatedAssemblyName() != null && !assNames.contains(e.getAssociatedAssemblyName())) ? e.getBreakpointQual() : 0);
			sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_SPLITREAD_QUAL, VcfFormatAttributes.BREAKPOINT_ASSEMBLED_NONSUPPORTING_SPLITREAD_QUAL,
					supportingSR, e -> (e.getAssociatedAssemblyName() != null && !assNames.contains(e.getAssociatedAssemblyName())) ? e.getBreakpointQual() : 0);
		}
		if (isUpdateAssemblyInformation()) {
			// Assembly breakdown
			int[] asr = new int[processContext.getCategoryCount()];
			int[] asrp = new int[processContext.getCategoryCount()];
			int[] basr = new int[processContext.getCategoryCount()];
			int[] basrp = new int[processContext.getCategoryCount()];
			for (int i = 0; i < processContext.getCategoryCount(); i++) {
				final int category = i;
				asr[category] = Stream.concat(Stream.concat(supportingAS.stream(), supportingRAS.stream()), supportingCAS.stream())
						.mapToInt(ass -> aaLookup.get(ass).getSupportingReadCount(
								ass.getBreakendAssemblyContigOffset(),
								ImmutableSet.of(category),
								ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read),
								(AssemblyEvidenceSource)ass.getEvidenceSource(), processContext))
						.sum();
				asrp[category] = Stream.concat(Stream.concat(supportingAS.stream(), supportingRAS.stream()), supportingCAS.stream())
						.mapToInt(ass -> aaLookup.get(ass).getSupportingReadCount(
								ass.getBreakendAssemblyContigOffset(),
								ImmutableSet.of(category),
								ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair),
								(AssemblyEvidenceSource)ass.getEvidenceSource(),
                                ass.source.getContext()))
						.sum();
				basr[category] = supportingBAS.stream()
						.mapToInt(ass -> aaLookup.get(ass).getSupportingReadCount(
								ass.getBreakendAssemblyContigOffset(),
								ImmutableSet.of(category),
								ImmutableSet.of(AssemblyEvidenceSupport.SupportType.Read),
								(AssemblyEvidenceSource)ass.getEvidenceSource(),
                                ass.source.getContext()))
						.sum();
				basrp[category] = supportingBAS.stream()
						.mapToInt(ass -> aaLookup.get(ass).getSupportingReadCount(
								ass.getBreakendAssemblyContigOffset(),
								ImmutableSet.of(category),
								ImmutableSet.of(AssemblyEvidenceSupport.SupportType.ReadPair),
								(AssemblyEvidenceSource)ass.getEvidenceSource(),
                                ass.source.getContext()))
						.sum();

				genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT.attribute(), asr[category]);
				genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT.attribute(), asrp[category]);
				genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_ASSEMBLY_READ_COUNT.attribute(), basr[category]);
				genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_ASSEMBLY_READPAIR_COUNT.attribute(), basrp[category]);
			}
			attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT.attribute(), IntStream.of(asr).sum());
			attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT.attribute(), IntStream.of(asrp).sum());
			attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_READ_COUNT.attribute(), IntStream.of(basr).sum());
			attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_READPAIR_COUNT.attribute(), IntStream.of(basrp).sum());
		}
		if (isUpdateAssemblyInformation() & isUpdateReadInformation()) {
			int[] supportingBreakpointFragments = new int[processContext.getCategoryCount()];
			int[] supportingBreakendFragments = new int[processContext.getCategoryCount()];
			for (int i = 0; i < processContext.getCategoryCount(); i++) {
				final int category = i;
				Set<String> bpfrags = supportingBreakpoint.stream()
						.map(e -> {
							return getOriginatingFragmentIDs(aaLookup, category, e);
						})
						.flatMap(x -> x.stream())
						.collect(Collectors.toSet());
				supportingBreakpointFragments[category] = bpfrags.size();
				Set<String> befrags = supportingBreakend.stream()
						.map(e -> {
							return getOriginatingFragmentIDs(aaLookup, category, e);
						})
						.flatMap(x -> x.stream())
						.collect(Collectors.toSet());
				befrags.removeAll(bpfrags);
				supportingBreakendFragments[category] = befrags.size();

				genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKPOINT_VARIANT_FRAGMENTS.attribute(), supportingBreakpointFragments[category]);
				genotypeBuilder.get(category).attribute(VcfFormatAttributes.BREAKEND_VARIANT_FRAGMENTS.attribute(), supportingBreakendFragments[category]);
			}
			attribute(VcfInfoAttributes.BREAKPOINT_VARIANT_FRAGMENTS.attribute(), IntStream.of(supportingBreakpointFragments).sum());
			attribute(VcfInfoAttributes.BREAKEND_VARIANT_FRAGMENTS.attribute(), IntStream.of(supportingBreakendFragments).sum());
		}
	}

	private void updateVariantQualityScoreAttributes(Map<SingleReadEvidence, AssemblyAttributes> aaLookup) {
		if (isUpdateAssemblyInformation() && isUpdateReadInformation()) {
			attribute(VcfInfoAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
			double beQual = supportingBAS.stream().mapToDouble(e -> e.getBreakendQual()).sum()
					+ supportingSC.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakendQual()).sum()
					+ supportingOEA.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakendQual()).sum();
			double bpQual = supportingAS.stream().mapToDouble(e -> ((DirectedBreakpoint) e).getBreakpointQual()).sum()
					+ supportingRAS.stream().mapToDouble(e -> ((DirectedBreakpoint) e).getBreakpointQual()).sum()
					+ supportingCAS.stream().mapToDouble(e -> ((DirectedBreakpoint) e).getBreakpointQual()).sum()
					+ supportingSR.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakpointQual()).sum()
					+ supportingIndel.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakpointQual()).sum()
					+ supportingDP.stream().flatMap(l -> l.stream()).mapToDouble(e -> e.getBreakpointQual()).sum();
			attribute(VcfInfoAttributes.BREAKEND_QUAL.attribute(), beQual);
			phredScore(isBreakend() ? beQual : bpQual);
		}
		double[] asq = null;
		double[] rasq = null;
		double[] casq = null;
		double[] basq = null;
		// Count
		if (isUpdateAssemblyInformation()) {
			attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT, supportingAS.size());
			attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE, supportingRAS.size());
			attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_COMPOUND, supportingCAS.size());
			attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_COUNT, supportingBAS.size());
			asq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_QUAL, supportingAS, aaLookup);
			rasq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, supportingRAS, aaLookup);
			casq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_COMPOUND, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_QUAL_COMPOUND, supportingCAS, aaLookup);
			basq = prorataAssemblyQualBreakdown(VcfInfoAttributes.BREAKEND_ASSEMBLY_QUAL, VcfFormatAttributes.BREAKEND_ASSEMBLY_QUAL, supportingBAS, aaLookup);
		}
		double[] srq = null;
		double[] iq = null;
		double[] rpq = null;
		double[] scq = null;
		double[] umq = null;
		if (isUpdateReadInformation()) {
			sumIntAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT, VcfFormatAttributes.BREAKPOINT_SPLITREAD_COUNT, supportingSR, e -> 1);
			sumIntAttr(VcfInfoAttributes.BREAKPOINT_INDEL_COUNT, VcfFormatAttributes.BREAKPOINT_INDEL_COUNT, supportingIndel, e -> 1);
			sumIntAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_READPAIR_COUNT, supportingDP, e -> 1);
			sumIntAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT, VcfFormatAttributes.BREAKEND_SOFTCLIP_COUNT, supportingSC, e -> 1);
			sumIntAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_COUNT, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_COUNT, supportingOEA, e -> 1);
			srq = sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_QUAL, VcfFormatAttributes.BREAKPOINT_SPLITREAD_QUAL, supportingSR, e -> e.getBreakpointQual());
			iq = sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_INDEL_QUAL, VcfFormatAttributes.BREAKPOINT_INDEL_QUAL, supportingIndel, e -> e.getBreakpointQual());
			rpq = sumDoubleAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_QUAL, VcfFormatAttributes.BREAKPOINT_READPAIR_QUAL, supportingDP, e -> e.getBreakpointQual());
			scq = sumDoubleAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_QUAL, VcfFormatAttributes.BREAKEND_SOFTCLIP_QUAL, supportingSC, e -> e.getBreakendQual());
			umq = sumDoubleAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_QUAL, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_QUAL, supportingOEA, e -> e.getBreakendQual());
			sumDoubleAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_QUAL, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_QUAL, supportingOEA, e -> e.getBreakendQual());
		}
		if (isUpdateAssemblyInformation() && isUpdateReadInformation()) {
			for (int i = 0; i < processContext.getCategoryCount(); i++) {
				genotypeBuilder.get(i).attribute(VcfFormatAttributes.BREAKPOINT_QUAL.attribute(), srq[i] + iq[i] + rpq[i] + asq[i] + rasq[i] + casq[i]);
				genotypeBuilder.get(i).attribute(VcfFormatAttributes.BREAKEND_QUAL.attribute(), scq[i] + umq[i] + basq[i]);
			}
		}
		if (isUpdateReadInformation()) {
			if (supportingBreakend.size() > 0) {
				attribute(VcfInfoAttributes.BREAKEND_MEAN_SUPPORTING_MAPQ, supportingBreakend.stream().mapToDouble(de -> de.getLocalMapq()).average().getAsDouble());
				attribute(VcfInfoAttributes.BREAKEND_MAX_SUPPORTING_MAPQ, supportingBreakend.stream().mapToDouble(de -> de.getLocalMapq()).max().getAsDouble());
				attribute(VcfInfoAttributes.BREAKEND_MIN_SUPPORTING_MAPQ, supportingBreakend.stream().mapToDouble(de -> de.getLocalMapq()).min().getAsDouble());
			}
			if (supportingBreakpoint.size() > 0) {
				attribute(VcfInfoAttributes.MEAN_SUPPORTING_MAPQ, supportingBreakpoint.stream().mapToDouble(de -> (de.getLocalMapq() + de.getRemoteMapq()) / 2).average().getAsDouble());
				attribute(VcfInfoAttributes.MAX_SUPPORTING_MAPQ, supportingBreakpoint.stream().mapToDouble(de -> Math.max(de.getLocalMapq(), de.getRemoteMapq())).max().getAsDouble());
				attribute(VcfInfoAttributes.MIN_SUPPORTING_MAPQ, supportingBreakpoint.stream().mapToDouble(de -> Math.min(de.getLocalMapq(), de.getRemoteMapq())).min().getAsDouble());
			}
		}
	}

	private Collection<String> getOriginatingFragmentIDs(Map<SingleReadEvidence, AssemblyAttributes> aaLookup, int category, DirectedEvidence e) {
		if (AssemblyAttributes.isAssembly(e)) {
			SingleReadEvidence ass = (SingleReadEvidence)e;
			return aaLookup.get(ass).getOriginatingFragmentID(Range.closed(ass.getBreakendAssemblyContigOffset(), ass.getBreakendAssemblyContigOffset()), ImmutableSet.of(category), null, (AssemblyEvidenceSource)e.getEvidenceSource());
		}
		return e.getOriginatingFragmentID(category);
	}

	/*
	 * To get per-sample assembly scores, we pro-rata the assembly score by the contribution from each category.
	 * 
	 * Note that the sum of the contributing breakend scores is not necessarily the same as the assembly breakpoint score.
	 * This means that each assembly must be individually pro-rataed. 
	 */
	private <T extends SingleReadEvidence> double[] prorataAssemblyQualBreakdown(
			VcfInfoAttributes infoAttr, VcfFormatAttributes attr, List<T> assemblies,
			Map<SingleReadEvidence, AssemblyAttributes> aaLookup) {
		double totalAssQual = 0;
		double[] prorata = new double[processContext.getCategoryCount()];
		for (SingleReadEvidence ass : assemblies) {
			int offset = ass.getBreakendAssemblyContigOffset();
			AssemblyAttributes aa = aaLookup.get(ass);
			double assQual = (ass instanceof DirectedBreakpoint) ? ((DirectedBreakpoint)ass).getBreakpointQual() : ass.getBreakendQual();
			double[] breakdownQual = new double[processContext.getCategoryCount()];
			for (int i = 0; i < breakdownQual.length; i++) {
				breakdownQual[i] = aa.getSupportingQualScore(offset, ImmutableSet.of(i), null, (AssemblyEvidenceSource)ass.getEvidenceSource(), processContext);
			}
			double breakdownTotal = DoubleStream.of(breakdownQual).sum();
			if (breakdownTotal != 0) { // defensive check to mitigate impact of 0 qual assemblies (#156)
				for (int i = 0; i < prorata.length; i++) {
					prorata[i] += assQual * (breakdownQual[i] / breakdownTotal);
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
		List<VcfFilter> filters = processContext.getVariantCallingParameters().calculateCommonFilters(variant);
		if (variant instanceof VariantContextDirectedBreakpoint) {
			filters.addAll(processContext.getVariantCallingParameters().calculateBreakpointFilters((VariantContextDirectedBreakpoint)variant));
		} else {
			filters.addAll(processContext.getVariantCallingParameters().calculateSingleBreakendFilters(variant));
		}
		if (!filters.isEmpty()) {
			VariantContextBuilder builder = new VariantContextBuilder(variant);
			for (VcfFilter f : filters) {
				builder.filter(f.filter());
			}
			variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext.getDictionary(), variant.source, builder.make());
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

	public boolean isUpdateAssemblyInformation() {
		return updateAssemblyInformation;
	}

	public void setUpdateAssemblyInformation(boolean updateAssemblyInformation) {
		this.updateAssemblyInformation = updateAssemblyInformation;
	}

	public boolean isUpdateReadInformation() {
		return updateReadInformation;
	}

	public void setUpdateReadInformation(boolean updateReadInformation) {
		this.updateReadInformation = updateReadInformation;
	}

	public boolean isUpdateVariantQualityScore() {
		return updateReadInformation && updateAssemblyInformation;
	}
}
