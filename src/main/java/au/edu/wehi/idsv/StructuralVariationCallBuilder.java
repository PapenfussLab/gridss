package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Lists;
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
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class StructuralVariationCallBuilder extends IdsvVariantContextBuilder {
	private static final Log log = Log.getInstance(StructuralVariationCallBuilder.class);
	private final ProcessingContext processContext;
	private final VariantContextDirectedEvidence parent;
	private final Set<String> encounteredEvidenceIDs;
	private final RangeSet<Integer> anchoredBases = TreeRangeSet.create();
	private final int categories;
	private DirectedBreakpoint bestExactBreakpoint = null;
	private int fBREAKPOINT_ASSEMBLY_COUNT;
	private int[] fBREAKPOINT_READPAIR_COUNT;
	private int[] fBREAKPOINT_SPLITREAD_COUNT;
	private int[] fBREAKPOINT_INDEL_COUNT;
	private int fBREAKPOINT_ASSEMBLY_COUNT_REMOTE;
	//private int[] fBREAKPOINT_SPLITREAD_COUNT_REMOTE;
	private int[] fBREAKPOINT_ASSEMBLY_READPAIR_COUNT;
	private int[] fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT;
	//private int[] fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT;
	//private int[] fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT;	
	private double fBREAKPOINT_ASSEMBLY_QUAL;
	private double[] fBREAKPOINT_READPAIR_QUAL;
	private double[] fBREAKPOINT_SPLITREAD_QUAL;
	private double[] fBREAKPOINT_INDEL_QUAL;
	private double fBREAKPOINT_ASSEMBLY_QUAL_REMOTE;
	//private double[] fBREAKPOINT_SPLITREAD_QUAL_REMOTE;
	private int fBREAKEND_ASSEMBLY_COUNT;
	private int[] fBREAKEND_UNMAPPEDMATE_COUNT;
	private int[] fBREAKEND_SOFTCLIP_COUNT;
	private double fBREAKEND_ASSEMBLY_QUAL;
	private double[] fBREAKEND_UNMAPPEDMATE_QUAL;
	private double[] fBREAKEND_SOFTCLIP_QUAL;
	private List<String> BREAKEND_ASSEMBLY_ID = Lists.newArrayList();
	private List<GenotypeBuilder> genotypeList;
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent) {
		this(processContext, parent, true);
	}
	public StructuralVariationCallBuilder(ProcessingContext processContext, VariantContextDirectedEvidence parent, boolean deduplicateEvidence) {
		super(processContext, parent);
		this.processContext = processContext;
		this.parent = parent;
		this.categories = processContext.getCategoryCount();
		this.encounteredEvidenceIDs = deduplicateEvidence ? new HashSet<String>() : null;
		fBREAKPOINT_ASSEMBLY_COUNT = 0;
		fBREAKPOINT_READPAIR_COUNT = new int[categories];
		fBREAKPOINT_SPLITREAD_COUNT = new int[categories];
		fBREAKPOINT_INDEL_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_COUNT_REMOTE = 0;
		//fBREAKPOINT_SPLITREAD_COUNT_REMOTE = new int[categories];
		fBREAKPOINT_ASSEMBLY_READPAIR_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT = new int[categories];
		//fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT = new int[categories];
		//fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_QUAL = 0;
		fBREAKPOINT_READPAIR_QUAL = new double[categories];
		fBREAKPOINT_SPLITREAD_QUAL = new double[categories];
		fBREAKPOINT_INDEL_QUAL = new double[categories];
		fBREAKPOINT_ASSEMBLY_QUAL_REMOTE = 0;
		//fBREAKPOINT_SPLITREAD_QUAL_REMOTE = new double[categories];
		fBREAKEND_ASSEMBLY_COUNT = 0;
		fBREAKEND_UNMAPPEDMATE_COUNT = new int[categories];
		fBREAKEND_SOFTCLIP_COUNT = new int[categories];
		fBREAKEND_ASSEMBLY_QUAL = 0;
		fBREAKEND_UNMAPPEDMATE_QUAL = new double[categories];
		fBREAKEND_SOFTCLIP_QUAL = new double[categories];
		genotypeList = IntStream.range(0, categories)
				.mapToObj(i -> new GenotypeBuilder(processContext.getCategoryLabel(i))
						.alleles(Arrays.asList(Allele.NO_CALL))
						.phased(false)
						.noAD()
						.noDP()
						.noGQ()
						.noPL())
				.collect(Collectors.toList());
	}
	private static int deduplicationMessageCount = 0;
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (!isSupportingEvidence(evidence)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s %s does not provide support for call at %s",
					evidence.getEvidenceID(),
					evidence.getBreakendSummary(),
					parent.getBreakendSummary()));
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
		if (evidence instanceof DirectedBreakpoint && evidence.isBreakendExact()) {
			DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
			if (ByBestBreakpointDesc.compare(bp, bestExactBreakpoint) < 0) {
				bestExactBreakpoint = (DirectedBreakpoint)evidence;
			}
		}
		int category = -1;
		if (evidence.getEvidenceSource() instanceof SAMEvidenceSource) {
			category = ((SAMEvidenceSource)evidence.getEvidenceSource()).getSourceCategory();
		}
		DirectedBreakpoint bp = (evidence instanceof DirectedBreakpoint) ? (DirectedBreakpoint)evidence : null;
		if (evidence instanceof NonReferenceReadPair && category >= 0) {
			if (evidence instanceof DiscordantReadPair) {
				fBREAKPOINT_READPAIR_COUNT[category]++;
				fBREAKPOINT_READPAIR_QUAL[category] += bp.getBreakpointQual();
				processAnchor(((NonReferenceReadPair)evidence).getLocalledMappedRead());
			} else {
				fBREAKEND_UNMAPPEDMATE_COUNT[category]++;
				fBREAKEND_UNMAPPEDMATE_QUAL[category] += evidence.getBreakendQual();
			}
		} else {
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			if (AssemblyAttributes.isAssembly(sre.getSAMRecord())) {
				AssemblyAttributes attr = new AssemblyAttributes(sre.getSAMRecord());
				if (sre instanceof DirectedBreakpoint) {
					for (int i = 0; i < categories; i++) {
						fBREAKPOINT_ASSEMBLY_READPAIR_COUNT[i] += attr.getAssemblySupportCountReadPair(i);
						fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT[i] += attr.getAssemblySupportCountSoftClip(i);
						//fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT[i] += attr.getAssemblyNonSupportingReadPairCount(i);
						//fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT[i]  += attr.getAssemblyNonSupportingSoftClipCount(i);
						// TODO quals, remotes, etc
					}
					if (sre.getSAMRecord().getSupplementaryAlignmentFlag() || attr.getAssemblyDirection() != sre.getBreakendSummary().direction) {
						// remote breakpoint assembly
						fBREAKPOINT_ASSEMBLY_COUNT_REMOTE++;
						fBREAKPOINT_ASSEMBLY_QUAL_REMOTE += bp.getBreakpointQual();
					} else {
						// local breakpoint assembly
						fBREAKPOINT_ASSEMBLY_COUNT++;
						fBREAKPOINT_ASSEMBLY_QUAL += bp.getBreakpointQual();
					}
					BREAKEND_ASSEMBLY_ID.add(sre.getSAMRecord().getReadName());
				} else {
					// breakend assembly
					fBREAKEND_ASSEMBLY_COUNT++;
					fBREAKEND_ASSEMBLY_QUAL += evidence.getBreakendQual();
				}
			} else if (category >= 0) {
				if (evidence instanceof SoftClipEvidence) {
					fBREAKEND_SOFTCLIP_COUNT[category]++;
					fBREAKEND_SOFTCLIP_QUAL[category] += evidence.getBreakendQual();
				} else if (evidence instanceof SplitReadEvidence) {
					fBREAKPOINT_SPLITREAD_COUNT[category]++;
					fBREAKPOINT_SPLITREAD_QUAL[category] += bp.getBreakpointQual();
				} else if (evidence instanceof IndelEvidence) {
					fBREAKPOINT_INDEL_COUNT[category]++;
					fBREAKPOINT_INDEL_QUAL[category] += bp.getBreakpointQual();
				}
			}
			if (sre instanceof DirectedBreakpoint) {
				processAnchor(sre.getSAMRecord());
			}
		}
		return this;
	}
	private void processAnchor(SAMRecord record) {
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
		return parent.getBreakendSummary().overlaps(bs);
	}
	private double sum(double[] v) {
		double total = 0;
		for (double x : v) {
			total += x;
		}
		return total;
	}
	private float[] tof(double[] v) {
		float[] f = new float[v.length];
		for (int i = 0; i < v.length; i++) {
			f[i] = (float)v[i];
		}
		return f;
	}
	private int sum(int[] v) {
		int total = 0;
		for (int x : v) {
			total += x;
		}
		return total;
	}
	public VariantContextDirectedEvidence make() {
		List<GenotypeBuilder> genotypeList = new ArrayList<>();
		for (int i = 0; i < processContext.getCategoryCount(); i++) {
			GenotypeBuilder gb = new GenotypeBuilder(processContext.getCategoryLabel(i))
					.alleles(Arrays.asList(Allele.NO_CALL))
					.phased(false)
					.noAD()
					.noDP()
					.noGQ()
					.noPL();
			//TODO:: add category attributes
			genotypeList.add(gb);
		}
		attribute(VcfInfoAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		attribute(VcfInfoAttributes.BREAKEND_QUAL.attribute(),
				fBREAKEND_ASSEMBLY_QUAL
				+ sum(fBREAKEND_SOFTCLIP_QUAL)
				+ sum(fBREAKEND_UNMAPPEDMATE_QUAL));
		phredScore(
				fBREAKPOINT_ASSEMBLY_QUAL
				+ fBREAKPOINT_ASSEMBLY_QUAL_REMOTE
				+ sum(fBREAKPOINT_READPAIR_QUAL)
				+ sum(fBREAKPOINT_SPLITREAD_QUAL)
				+ sum(fBREAKPOINT_INDEL_QUAL));
				//+ sum(fBREAKPOINT_SPLITREAD_QUAL_REMOTE));
		
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT, fBREAKPOINT_ASSEMBLY_COUNT);
		sumAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_READPAIR_COUNT, fBREAKPOINT_READPAIR_COUNT);
		sumAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT, VcfFormatAttributes.BREAKPOINT_SPLITREAD_COUNT, fBREAKPOINT_SPLITREAD_COUNT);
		sumAttr(VcfInfoAttributes.BREAKPOINT_INDEL_COUNT, VcfFormatAttributes.BREAKPOINT_INDEL_COUNT, fBREAKPOINT_INDEL_COUNT);
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE, fBREAKPOINT_ASSEMBLY_COUNT_REMOTE);
		sumAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, fBREAKPOINT_ASSEMBLY_READPAIR_COUNT);
		sumAttr(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT, VcfFormatAttributes.BREAKPOINT_ASSEMBLY_READ_COUNT, fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT);
		//sumAttr(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT, fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT);
		//sumAttr(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READ_COUNT, fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT);
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL, (float)fBREAKPOINT_ASSEMBLY_QUAL);
		sumAttr(VcfInfoAttributes.BREAKPOINT_READPAIR_QUAL, VcfFormatAttributes.BREAKPOINT_READPAIR_QUAL, fBREAKPOINT_READPAIR_QUAL);
		sumAttr(VcfInfoAttributes.BREAKPOINT_SPLITREAD_QUAL, VcfFormatAttributes.BREAKPOINT_SPLITREAD_QUAL, fBREAKPOINT_SPLITREAD_QUAL);
		sumAttr(VcfInfoAttributes.BREAKPOINT_INDEL_QUAL, VcfFormatAttributes.BREAKPOINT_INDEL_QUAL, fBREAKPOINT_INDEL_QUAL);
		attribute(VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, (float)fBREAKPOINT_ASSEMBLY_QUAL_REMOTE);
		attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_COUNT, fBREAKEND_ASSEMBLY_COUNT);
		sumAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_COUNT, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_COUNT, fBREAKEND_UNMAPPEDMATE_COUNT);
		sumAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT, VcfFormatAttributes.BREAKEND_SOFTCLIP_COUNT, fBREAKEND_SOFTCLIP_COUNT);
		attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_QUAL, (float)fBREAKEND_ASSEMBLY_QUAL);
		sumAttr(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_QUAL, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_QUAL, fBREAKEND_UNMAPPEDMATE_QUAL);
		sumAttr(VcfInfoAttributes.BREAKEND_SOFTCLIP_QUAL, VcfFormatAttributes.BREAKEND_SOFTCLIP_QUAL, fBREAKEND_SOFTCLIP_QUAL);
		attribute(VcfInfoAttributes.BREAKEND_ASSEMBLY_ID, BREAKEND_ASSEMBLY_ID.toArray(new String[0]));
		
		genotypes(genotypeList.stream().map(gb -> gb.make()).collect(Collectors.toList()));
			
		
		String untemplated = parent.getBreakpointSequenceString();
		String homo = "";
		if (bestExactBreakpoint != null) {
			untemplated = bestExactBreakpoint.getUntemplatedSequence();
			breakpoint(bestExactBreakpoint.getBreakendSummary(), untemplated);
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
		attribute(VcfInfoAttributes.SUPPORT_CIGAR, makeCigar(anchoredBases, bestExactBreakpoint != null ? bestExactBreakpoint.getBreakendSummary() : parent.getBreakendSummary()).toString());
		// id(parent.getID()); // can't change from parent ID as the id is already referenced in the MATEID of the other breakend  
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = applyFilters(variant);
		//variant = Models.calculateSomatic(variant);
		//assert(sanitycheck(variant));
		return variant;
	}
	private void sumAttr(VcfInfoAttributes infoAttr, VcfFormatAttributes formatAttr, int[] values) {
		attribute(infoAttr, sum(values));
		for (int i = 0; i < values.length; i++) {
			genotypeList.get(i).attribute(formatAttr.name(), values[i]);
		}
	}
	private void sumAttr(VcfInfoAttributes infoAttr, VcfFormatAttributes formatAttr, double[] values) {
		attribute(infoAttr, (float)sum(values));
		for (int i = 0; i < values.length; i++) {
			genotypeList.get(i).attribute(formatAttr.name(), (float)values[i]);
		}
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
