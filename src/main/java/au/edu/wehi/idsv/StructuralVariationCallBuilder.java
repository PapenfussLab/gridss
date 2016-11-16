package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import au.edu.wehi.idsv.alignment.BreakpointHomology;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;
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
	private final RangeSet<Integer> anchoredBases = TreeRangeSet.create();
	private final int categories;
	private DirectedBreakpoint bestExactBreakpoint = null;
	private int fBREAKPOINT_ASSEMBLY_COUNT;
	private int[] fBREAKPOINT_READPAIR_COUNT;
	private int[] fBREAKPOINT_SPLITREAD_COUNT;
	private int fBREAKPOINT_ASSEMBLY_COUNT_REMOTE;
	private int[] fBREAKPOINT_SPLITREAD_COUNT_REMOTE;
	private int[] fBREAKPOINT_ASSEMBLY_READPAIR_COUNT;
	private int[] fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT;
	private int[] fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT;
	private int[] fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT;	
	private double fBREAKPOINT_ASSEMBLY_QUAL;
	private double[] fBREAKPOINT_READPAIR_QUAL;
	private double[] fBREAKPOINT_SPLITREAD_QUAL;
	private double fBREAKPOINT_ASSEMBLY_QUAL_REMOTE;
	private double[] fBREAKPOINT_SPLITREAD_QUAL_REMOTE;
	private int fBREAKEND_ASSEMBLY_COUNT;
	private int[] fBREAKEND_UNMAPPEDMATE_COUNT;
	private int[] fBREAKEND_SOFTCLIP_COUNT;
	private double fBREAKEND_ASSEMBLY_QUAL;
	private double[] fBREAKEND_UNMAPPEDMATE_QUAL;
	private double[] fBREAKEND_SOFTCLIP_QUAL;
	private List<String> BREAKEND_ASSEMBLY_ID = Lists.newArrayList();

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
		fBREAKPOINT_ASSEMBLY_COUNT_REMOTE = 0;
		fBREAKPOINT_SPLITREAD_COUNT_REMOTE = new int[categories];
		fBREAKPOINT_ASSEMBLY_READPAIR_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT = new int[categories];
		fBREAKPOINT_ASSEMBLY_QUAL = 0;
		fBREAKPOINT_READPAIR_QUAL = new double[categories];
		fBREAKPOINT_SPLITREAD_QUAL = new double[categories];
		fBREAKPOINT_ASSEMBLY_QUAL_REMOTE = 0;
		fBREAKPOINT_SPLITREAD_QUAL_REMOTE = new double[categories];
		fBREAKEND_ASSEMBLY_COUNT = 0;
		fBREAKEND_UNMAPPEDMATE_COUNT = new int[categories];
		fBREAKEND_SOFTCLIP_COUNT = new int[categories];
		fBREAKEND_ASSEMBLY_QUAL = 0;
		fBREAKEND_UNMAPPEDMATE_QUAL = new double[categories];
		fBREAKEND_SOFTCLIP_QUAL = new double[categories];
	}
	public StructuralVariationCallBuilder addEvidence(DirectedEvidence evidence) {
		if (evidence == null) throw new NullPointerException();
		if (!isSupportingEvidence(evidence)) {
			throw new IllegalArgumentException(String.format("Sanity check failure: Evidence %s does not provide support for call at %s", evidence.getBreakendSummary(), parent.getBreakendSummary()));
		}
		String eid = evidence.getEvidenceID();
		if (encounteredEvidenceIDs != null) {
			if (encounteredEvidenceIDs.contains(eid)) {
				log.debug(String.format("Deduplicating %s from %s", eid, parent.getID()));
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
		int category = evidence.getEvidenceSource() instanceof SAMEvidenceSource ? ((SAMEvidenceSource)evidence.getEvidenceSource()).getSourceCategory() : -1;
		if (evidence instanceof RealignedRemoteSAMRecordAssemblyEvidence) {
			add((RealignedRemoteSAMRecordAssemblyEvidence)evidence);
		} else if (evidence instanceof SpanningSAMRecordAssemblyEvidence) {
			add((SpanningSAMRecordAssemblyEvidence)evidence);
		} else if (evidence instanceof RealignedSAMRecordAssemblyEvidence) {
			add((RealignedSAMRecordAssemblyEvidence)evidence);
		} else if (evidence instanceof SAMRecordAssemblyEvidence) {
			add((SAMRecordAssemblyEvidence)evidence);
		} else if (evidence instanceof DiscordantReadPair) {
			add((DiscordantReadPair)evidence, category);
		} else if (evidence instanceof UnmappedMateReadPair) {
			add((UnmappedMateReadPair)evidence, category);
		}  else if (evidence instanceof RealignedRemoteSoftClipEvidence) {
			add((RealignedRemoteSoftClipEvidence)evidence, category);
		} else if (evidence instanceof RealignedSoftClipEvidence) {
			add((RealignedSoftClipEvidence)evidence, category);
		} else if (evidence instanceof SoftClipEvidence) {
			add((SoftClipEvidence)evidence, category);
		} else {
			throw new IllegalArgumentException(String.format("Unknown evidence type %s", evidence.getClass()));
		}
		return this;
	}
	private void addAssemblyCategories(SAMRecordAssemblyEvidence evidence) {
		for (int i = 0; i < categories; i++) {
			fBREAKPOINT_ASSEMBLY_READPAIR_COUNT[i] += evidence.getAssemblySupportCountReadPair(i);
			fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT[i] += evidence.getAssemblySupportCountSoftClip(i);
			fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT[i] += evidence.getAssemblyNonSupportingReadPairCount(i);
			fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT[i]  += evidence.getAssemblyNonSupportingSoftClipCount(i);
			// TODO quals, remotes, etc
		}
	}
	private void add(RealignedRemoteSAMRecordAssemblyEvidence evidence) {
		fBREAKPOINT_ASSEMBLY_COUNT_REMOTE++;
		fBREAKPOINT_ASSEMBLY_QUAL_REMOTE += evidence.getBreakpointQual();
		BREAKEND_ASSEMBLY_ID.add(evidence.asLocal().getEvidenceID());
		addAssemblyCategories(evidence);
		if (evidence.isBreakendExact()) {
			processAnchor(evidence.getRemoteSAMRecord());
		}
	}
	private void add(SpanningSAMRecordAssemblyEvidence evidence) {
		BreakendDirection direction = evidence.getParentAssembly().getBreakendDirection();
		if (direction != evidence.getBreakendSummary().direction) {
			fBREAKPOINT_ASSEMBLY_COUNT_REMOTE++;
			fBREAKPOINT_ASSEMBLY_QUAL_REMOTE += evidence.getBreakpointQual();
		}
		if (direction != evidence.getBreakendSummary().direction.reverse()) {
			fBREAKPOINT_ASSEMBLY_COUNT++;
			fBREAKPOINT_ASSEMBLY_QUAL += evidence.getBreakpointQual();
		}
		BREAKEND_ASSEMBLY_ID.add(evidence.getEvidenceID());
		addAssemblyCategories(evidence);
		processAnchor(evidence.getSAMRecord());
	}
	private void add(RealignedSAMRecordAssemblyEvidence evidence) {
		fBREAKPOINT_ASSEMBLY_COUNT++;
		fBREAKPOINT_ASSEMBLY_QUAL += evidence.getBreakpointQual();
		BREAKEND_ASSEMBLY_ID.add(evidence.getEvidenceID());
		addAssemblyCategories(evidence);
		processAnchor(evidence.getSAMRecord());
	}
	private void add(SAMRecordAssemblyEvidence evidence) {
		fBREAKEND_ASSEMBLY_COUNT++;
		fBREAKEND_ASSEMBLY_QUAL += evidence.getBreakendQual();
		processAnchor(evidence.getSAMRecord());
	}
	private void add(DiscordantReadPair evidence, int category) {
		fBREAKPOINT_READPAIR_COUNT[category]++;
		fBREAKPOINT_READPAIR_QUAL[category] += evidence.getBreakpointQual();
		processAnchor(evidence.getLocalledMappedRead());
	}
	private void add(UnmappedMateReadPair evidence, int category) {
		fBREAKEND_UNMAPPEDMATE_COUNT[category]++;
		fBREAKEND_UNMAPPEDMATE_QUAL[category] += evidence.getBreakendQual();
		processAnchor(evidence.getLocalledMappedRead());
	}
	private void add(RealignedRemoteSoftClipEvidence evidence, int category) {
		fBREAKPOINT_SPLITREAD_COUNT_REMOTE[category]++;
		fBREAKPOINT_SPLITREAD_QUAL_REMOTE[category] += evidence.getBreakpointQual();
		processAnchor(evidence.getSAMRecord());
	}
	private void add(RealignedSoftClipEvidence evidence, int category) {
		fBREAKPOINT_SPLITREAD_COUNT[category]++;
		fBREAKPOINT_SPLITREAD_QUAL[category] += evidence.getBreakpointQual();
		processAnchor(evidence.getSAMRecord());
	}
	private void add(SoftClipEvidence evidence, int category) {
		fBREAKEND_SOFTCLIP_COUNT[category]++;
		fBREAKEND_SOFTCLIP_QUAL[category] += evidence.getBreakendQual();
		processAnchor(evidence.getSAMRecord());
	}
	private void processAnchor(SAMRecord record) {
		if (record.getReadUnmappedFlag()) return;
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
	public VariantContextDirectedEvidence make() {
		//extractAssemblySupport(); // don't extract - as extraction inflates evidence beyond original call
		attribute(VcfAttributes.CALLED_QUAL.attribute(), parent.getPhredScaledQual());
		attribute(VcfAttributes.BREAKEND_QUAL.attribute(),
				fBREAKEND_ASSEMBLY_QUAL
				+ sum(fBREAKEND_SOFTCLIP_QUAL)
				+ sum(fBREAKEND_UNMAPPEDMATE_QUAL));
		phredScore(
				fBREAKPOINT_ASSEMBLY_QUAL
				+ fBREAKPOINT_ASSEMBLY_QUAL_REMOTE
				+ sum(fBREAKPOINT_READPAIR_QUAL)
				+ sum(fBREAKPOINT_SPLITREAD_QUAL)
				+ sum(fBREAKPOINT_SPLITREAD_QUAL_REMOTE));
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT, fBREAKPOINT_ASSEMBLY_COUNT);
		attribute(VcfAttributes.BREAKPOINT_READPAIR_COUNT, fBREAKPOINT_READPAIR_COUNT);
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT, fBREAKPOINT_SPLITREAD_COUNT);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE, fBREAKPOINT_ASSEMBLY_COUNT_REMOTE);
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_COUNT_REMOTE, fBREAKPOINT_SPLITREAD_COUNT_REMOTE);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_READPAIR_COUNT, fBREAKPOINT_ASSEMBLY_READPAIR_COUNT);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_SPLITREAD_COUNT, fBREAKPOINT_ASSEMBLY_SPLITREAD_COUNT);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT, fBREAKPOINT_ASSEMBLY_CONSCRIPTED_READPAIR_COUNT);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT, fBREAKPOINT_ASSEMBLY_CONSCRIPTED_SPLITREAD_COUNT);
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL, (float)fBREAKPOINT_ASSEMBLY_QUAL);
		attribute(VcfAttributes.BREAKPOINT_READPAIR_QUAL, tof(fBREAKPOINT_READPAIR_QUAL));
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL, tof(fBREAKPOINT_SPLITREAD_QUAL));
		attribute(VcfAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE, (float)fBREAKPOINT_ASSEMBLY_QUAL_REMOTE);
		attribute(VcfAttributes.BREAKPOINT_SPLITREAD_QUAL_REMOTE, tof(fBREAKPOINT_SPLITREAD_QUAL_REMOTE));
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_COUNT, fBREAKEND_ASSEMBLY_COUNT);
		attribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_COUNT, fBREAKEND_UNMAPPEDMATE_COUNT);
		attribute(VcfAttributes.BREAKEND_SOFTCLIP_COUNT, fBREAKEND_SOFTCLIP_COUNT);
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_QUAL, (float)fBREAKEND_ASSEMBLY_QUAL);
		attribute(VcfAttributes.BREAKEND_UNMAPPEDMATE_QUAL, tof(fBREAKEND_UNMAPPEDMATE_QUAL));
		attribute(VcfAttributes.BREAKEND_SOFTCLIP_QUAL, tof(fBREAKEND_SOFTCLIP_QUAL));
		attribute(VcfAttributes.BREAKEND_ASSEMBLY_ID, BREAKEND_ASSEMBLY_ID.toArray(new String[0]));
		
		BreakendSummary bs = parent.getBreakendSummary();
		String untemplated = parent.getBreakpointSequenceString();
		String homo = "";
		if (bestExactBreakpoint != null) {
			bs = bestExactBreakpoint.getBreakendSummary();
			untemplated = bestExactBreakpoint.getUntemplatedSequence();
			breakpoint(bestExactBreakpoint.getBreakendSummary(), untemplated);
			homo = bestExactBreakpoint.getHomologySequence();
			rmAttribute(VcfSvConstants.IMPRECISE_KEY);
			
			BreakpointHomology bh = BreakpointHomology.calculate(
					processContext.getReference(),
					bestExactBreakpoint.getBreakendSummary(),
					untemplated,
					processContext.getVariantCallingParameters().maxBreakendHomologyLength,
					processContext.getVariantCallingParameters().breakendHomologyAlignmentMargin);
			int[] bounds;
			if (bs.direction == BreakendDirection.Forward) {
				bounds = new int[] { -bh.getLocalHomologyLength(), bh.getRemoteHomologyLength() };
			} else {
				bounds = new int[] { -bh.getRemoteHomologyLength(), bh.getLocalHomologyLength() };
			}
			attribute(VcfAttributes.INEXACT_HOMPOS, bounds);
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
		attribute(VcfAttributes.SUPPORT_CIGAR, makeCigar(anchoredBases, bestExactBreakpoint != null ? bestExactBreakpoint.getBreakendSummary() : parent.getBreakendSummary()).toString());
		// id(parent.getID()); // can't change from parent ID as the id is already referenced in the MATEID of the other breakend  
		VariantContextDirectedEvidence variant = (VariantContextDirectedEvidence)IdsvVariantContext.create(processContext, null, super.make());
		variant = applyFilters(variant);
		//variant = Models.calculateSomatic(variant);
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
					.compare(right.getBreakpointQual(), left.getBreakpointQual()) // desc
					// Take the one that aligned the most bases (more likely to remove REF FP calls)
					.compare(left.getUntemplatedSequence().length(), right.getUntemplatedSequence().length())
					// Ensure both sides of the breakpoint choose the same evidence in case of a tie 
					.compare(asLocal(left).getEvidenceID(), asLocal(right).getEvidenceID())
					.result();
		}
	}.nullsLast();
	private static DirectedBreakpoint asLocal(DirectedBreakpoint bp) {
		if (bp instanceof RemoteEvidence) {
			return ((RemoteEvidence) bp).asLocal();
		}
		return bp;
	}
}
