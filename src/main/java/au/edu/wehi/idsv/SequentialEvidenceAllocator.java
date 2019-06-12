package au.edu.wehi.idsv;

import au.edu.wehi.idsv.util.RangeUtil;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import au.edu.wehi.idsv.visualisation.TrackedBuffer;
import com.google.common.collect.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.lang3.StringUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Adds evidence supporting to each variant call. Both the variant calls and
 * evidence are required to be in-order.
 * 
 * @author Daniel Cameron
 *
 */
public class SequentialEvidenceAllocator implements Iterator<SequentialEvidenceAllocator.VariantEvidenceSupport>, TrackedBuffer {
	private static final Log log = Log.getInstance(SequentialEvidenceAllocator.class);
	private ProgressLogger progressLogger = new ProgressLogger(log);
	private final ProcessingContext context;
	private final int maxCallRange;
	private final boolean assignEvidenceToSingleBreakpoint;
	private final boolean preferToAssignBreakendReadsToVariantContainingAssembly = true;
	private final PeekingIterator<? extends DirectedEvidence> readIt;
	private final PeekingIterator<? extends DirectedEvidence> assemblyIt;
	private final Iterator<? extends VariantContextDirectedEvidence> callIt;
	private final OverlapLookup breakpointLookup;
	private final OverlapLookup breakendLookup;
	private final ArrayDeque<VariantEvidenceSupport> variantBuffer = new ArrayDeque<VariantEvidenceSupport>();
	private final Map<String, VariantEvidenceSupport> bufferedVariantId = new HashMap<String, VariantEvidenceSupport>();
	private final SetMultimap<String, VariantEvidenceSupport> assemblyAllocationLookup = HashMultimap.create();
	public class VariantEvidenceSupport {
		private final String id;
		private final String parid;
		private final String eventid;
		private final long startLocation;
		//public final long endLocation;
		private final BreakendSummary location;
		private final float score;
		public final List<DirectedEvidence> support = new ArrayList<>();
		public final VariantContextDirectedEvidence variant;
		private VariantEvidenceSupport(VariantContextDirectedEvidence call) {
			this.variant = call;
			this.id = call.hasID() ? call.getID() : null;
			this.parid = call.hasID() ? (String)call.getAttribute(VcfSvConstants.PARTNER_BREAKEND_ID_KEY, null) : null;
			this.eventid = (String)call.getAttribute(VcfSvConstants.BREAKEND_EVENT_ID_KEY, null);
			this.score = (float)call.getPhredScaledQual();
			assert(this.score >= 0); // variant must have score set
			this.location = call.getBreakendSummary();
			this.startLocation = context.getLinear().getStartLinearCoordinate(this.location);
			//this.endLocation = context.getLinear().getEndLinearCoordinate(this.location);
		}
		private void attributeEvidence(DirectedEvidence e) {
			support.add(e);
		}
		public String toString() {
			return String.format("%s %f %s", location, score, id);
		}
		@Override
		public int hashCode() {
			return id.hashCode();
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (!(obj instanceof VariantEvidenceSupport)) return false;
			VariantEvidenceSupport ves = (VariantEvidenceSupport)obj;
			if (id != null) return id.equals(ves.id);
			return location.equals(ves.location) && variant.equals(ves.variant);
		}
	}
	/**
	 * Orders variants by their score then position
	 * Positional comparison that returns the same order for both high and low breakends
	 * is required to ensure both sides of paired evidence is assigned to corresponding
	 * breakends of the same event. 
	 */
	private static final Ordering<VariantEvidenceSupport> ByScoreAscPositionDesc = new Ordering<VariantEvidenceSupport>() {
		public int compare(VariantEvidenceSupport o1, VariantEvidenceSupport o2) {
			ComparisonChain chain = ComparisonChain.start()
			        .compare(o1.score, o2.score);
			if (o1.location instanceof BreakpointSummary && o2.location instanceof BreakpointSummary) {
				chain = chain.compare((BreakpointSummary)o2.location, (BreakpointSummary)o1.location, BreakpointSummary.ByLowHigh);
			} else {
				chain = chain.compare(o2.location, o1.location, BreakendSummary.ByStartEnd);
			}
			chain = chain
			        .compare(o1.eventid, o2.eventid)
			        .compare(o1.id, o2.id);
			return chain.result();
		}
	};
	/**
	 * Creates an evidence allocator
	 * @param context processing context
	 * @param calls variant calls ordered by position
	 * @param evidence evidence order by breakend position
	 * @param maxCallWindowSize
	 * @param assignEvidenceToSingleBreakpoint uniquely assign evidence to only the highest scoring call
	 */
	public SequentialEvidenceAllocator(
			ProcessingContext context,
			Iterator<? extends VariantContextDirectedEvidence> calls,
			Iterator<? extends DirectedEvidence> reads,
			Iterator<? extends DirectedEvidence> assemblies,
			int maxCallWindowSize,
			boolean assignEvidenceToSingleBreakpoint) {
		this.context = context;
		this.maxCallRange = maxCallWindowSize;
		this.callIt = calls;
		this.readIt = Iterators.peekingIterator(reads);
		this.assemblyIt = Iterators.peekingIterator(assemblies);
		this.assignEvidenceToSingleBreakpoint = assignEvidenceToSingleBreakpoint;
		if (assignEvidenceToSingleBreakpoint) {
			// RemoteOverlapLookup only tracks the best local breakend so we can't use it if we want to assign to all matching breakpoints
			this.breakpointLookup = new RemoteOverlapLookup(this.context.getDictionary().getSequences().size());
		} else {
			this.breakpointLookup = new LocalOverlapLookup(this.context.getDictionary().getSequences().size());
		}
		this.breakendLookup = new LocalOverlapLookup(this.context.getDictionary().getSequences().size());
	}
	private void buffer(VariantContextDirectedEvidence variant) {
		VariantEvidenceSupport av = new VariantEvidenceSupport(variant);
		variantBuffer.add(av);
		if (StringUtils.isNotBlank(av.id)) {
			bufferedVariantId.put(av.id, av);
		}
		if (variant instanceof VariantContextDirectedBreakpoint) {
			breakpointLookup.add(av);
		} else {
			assert(!(variant.getBreakendSummary() instanceof BreakpointSummary)); // sanity check
			breakendLookup.add(av);
		}
	}
	@Override
	public boolean hasNext() {
		boolean hasNext = !variantBuffer.isEmpty() || callIt.hasNext();
		if (!hasNext) {
			if (Defaults.SANITY_CHECK_ITERATORS) {
				// we have no more calls to make so this doesn't actually need to be done
				// unless we're sanity checking
				while (assemblyIt.hasNext()) {
					assignEvidence(assemblyIt.next());
				}
				while (readIt.hasNext()) {
					assignEvidence(readIt.next());
				}
			}
		}
		return hasNext;
	}
	@Override
	public VariantEvidenceSupport next() {
		if (!hasNext()) throw new NoSuchElementException();
		if (variantBuffer.isEmpty()) {
			buffer(callIt.next());
		}
		VariantEvidenceSupport variant = variantBuffer.peek();
		// load all variant calls that could overlap an assembly
		bufferVariantsBefore(variant.startLocation + 3 * (maxCallRange + 1));
		// load all assemblies that could overlap the reads
		processEvidenceBefore(assemblyIt, variant.startLocation + 2 * (maxCallRange + 1));
		// load all reads that could overlap the variant call we're about to emit
		processEvidenceBefore(readIt, variant.startLocation + maxCallRange + 1);
		// Why do the above? 
		// We want to preferentially assign breakend reads to the variant
		// that they got assembled into. This means we need to load our
		// assemblies before our reads
		variant = variantBuffer.poll();
		if (StringUtils.isNotBlank(variant.id)) {
			bufferedVariantId.remove(variant.id);
		}
		if (variant.location instanceof BreakpointSummary) {
			breakpointLookup.remove(variant);
		} else {
			breakendLookup.remove(variant);
		}
		for (DirectedEvidence ass : variant.support) {
			if (AssemblyAttributes.isAssembly(ass)) {
				if (!assemblyAllocationLookup.remove(ass.getAssociatedAssemblyName(), variant) && assignEvidenceToSingleBreakpoint) {
					log.debug("Sanity failure: failed to remove assembly from lookup. Multiple evidence from single assembly assigned to this variant?");
				}
			}
		}
		return variant;
	}
	private void processEvidenceBefore(PeekingIterator<? extends DirectedEvidence> it, long position) {
		while (it.hasNext() && context.getLinear().getStartLinearCoordinate(it.peek().getBreakendSummary()) - context.getVariantCallingParameters().breakendMargin <= position) {
			assignEvidence(it.next());
		}
	}
	/**
	 * Assigns the given evidence to the appropriate call
	 * 
	 */
	private void assignEvidence(DirectedEvidence evidence) {
		BreakendSummary bs = evidence.getBreakendSummary();
		bs = context.getVariantCallingParameters().withMargin(bs);
		if (assignEvidenceToSingleBreakpoint) {
			assignToBest(bs, evidence);
		} else {
			assignToAll(bs, evidence);
		}
		if (evidence instanceof NonReferenceReadPair) {
			progressLogger.record(((NonReferenceReadPair)evidence).getLocalledMappedRead());
		} else if (evidence instanceof SingleReadEvidence) {
			progressLogger.record(((SingleReadEvidence)evidence).getSAMRecord());
		}
	}
	private void assignToBest(BreakendSummary bs, DirectedEvidence evidence) {
		VariantEvidenceSupport assignedTo;
		if (AssemblyAttributes.isAssembly(evidence)) {
			// breakend assemblies never get assigned to a breakpoint
			if (evidence instanceof DirectedBreakpoint) {
				assignedTo = assignToBestBreakpoint(bs, evidence);
			} else {
				assignedTo = assignToBestBreakend(bs, evidence);
			}
			if (assignedTo != null) {
				assemblyAllocationLookup.put(evidence.getAssociatedAssemblyName(), assignedTo);
			}
		} else {
			if (evidence instanceof DirectedBreakpoint) {
				// breakpoint evidence goes to the best breakpoint
				assignToBestBreakpoint(bs, evidence);
				// TODO should we attempt to fall back to breakend assignment?
			} else {
				// breakend evidence follows the assembly (if possible)
				VariantEvidenceSupport bestAssTo = null;
				if (preferToAssignBreakendReadsToVariantContainingAssembly && evidence.getAssociatedAssemblyName() != null) {
					Collection<VariantEvidenceSupport> hits = assemblyAllocationLookup.get(evidence.getAssociatedAssemblyName());
					if (hits != null) {
						for (VariantEvidenceSupport ves : hits) {
							if (ves.location.overlaps(bs) && (bestAssTo == null || bestAssTo.score < ves.score)) {
								bestAssTo = ves;
							}
						}
					}
				}
				if (bestAssTo != null) {
					bestAssTo.attributeEvidence(evidence);
				} else {
					if (assignToBestBreakpoint(bs, evidence) == null) {
						assignToBestBreakend(bs, evidence);
					}
				}
			}
		}
	}
	private VariantEvidenceSupport assignToBestBreakpoint(BreakendSummary bs, DirectedEvidence evidence) {
		VariantEvidenceSupport best = breakpointLookup.findBestOverlapping(bs);
		if (best != null) {
			VariantEvidenceSupport mate = bufferedVariantId.get(best.parid);
			if (mate != null && mate.location.overlaps(bs) && allocateToHighBreakend(evidence)) {
				// special case: matches both sides of the breakpoint 
				mate.attributeEvidence(evidence);
			} else {
				best.attributeEvidence(evidence);
			}
		}
		return best;
	}
	private VariantEvidenceSupport assignToBestBreakend(BreakendSummary bs, DirectedEvidence evidence) {
		VariantEvidenceSupport best = breakendLookup.findBestOverlapping(bs);
		if (best != null) {
			best.attributeEvidence(evidence);
		}
		return best;
	}
	private void assignToAll(BreakendSummary bs, DirectedEvidence evidence) {
		Iterator<VariantEvidenceSupport> it = breakpointLookup.findAllOverlapping(bs);
		if (!(evidence instanceof DirectedBreakpoint)) {
			it = Iterators.concat(it, breakendLookup.findAllOverlapping(bs));
		}
		while (it.hasNext()) {
			VariantEvidenceSupport v = it.next();
			assert(v.location.overlaps(bs));
			v.attributeEvidence(evidence);
		}
	}
	/**
	 * Determines which breakend to allocate evidence that overlaps both sides of the breakend
	 * @param e
	 * @return
	 */
	private boolean allocateToHighBreakend(DirectedEvidence evidence) {
		// want even allocation of evidence to both sides
		// want local and remote for same evidence allocated to different sides
		// want read pair evidence to be allocated to different sides
		String commonIdentifier;
		boolean flip = false;
		if (evidence instanceof NonReferenceReadPair) {
			// read name same for both sides of Discordant pairs
			SAMRecord local = ((NonReferenceReadPair) evidence).getLocalledMappedRead();
			commonIdentifier = local.getReadName();
			flip = local.getSecondOfPairFlag();
		} else if (evidence instanceof SingleReadEvidence) {
			SAMRecord local = ((SingleReadEvidence)evidence).getSAMRecord();
			commonIdentifier = local.getReadName();
			flip = local.getSupplementaryAlignmentFlag();
			if (evidence instanceof IndelEvidence) {
				flip ^= ((IndelEvidence)evidence).getBreakendSummary().isHighBreakend();
			}
		} else if (evidence instanceof VariantContextDirectedEvidence) {
			commonIdentifier = ((VariantContextDirectedEvidence)evidence).getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, evidence.getEvidenceID());
		} else {
			commonIdentifier = evidence.getEvidenceID();
		}
		boolean allocateLow = (Integer.bitCount(commonIdentifier.hashCode()) & 1) == 1; // randomly allocate high/low based on string hash
		allocateLow ^= flip;
		return allocateLow;
	}
	private void bufferVariantsBefore(long position) {
		while (callIt.hasNext() && (variantBuffer.isEmpty() || variantBuffer.peekLast().startLocation <= position)) {
			buffer(callIt.next());
		}
	}
	private String trackedBufferName_variantBuffer = "annotate.variantBuffer";
	private String trackedBufferName_bufferedVariantId = "annotate.bufferedVariantId";
	@Override
	public void setTrackedBufferContext(String context) {
		this.trackedBufferName_variantBuffer = context + ".annotate.variantBuffer";
		this.trackedBufferName_bufferedVariantId = context + ".annotate.bufferedVariantId";
	}
	@Override
	public List<NamedTrackedBuffer> currentTrackedBufferSizes() {
		return ImmutableList.of(
				new NamedTrackedBuffer(trackedBufferName_variantBuffer, variantBuffer.size()),
				new NamedTrackedBuffer(trackedBufferName_bufferedVariantId, bufferedVariantId.size())
				);
	}
	private static abstract class OverlapLookup {
		public abstract void add(VariantEvidenceSupport ves);
		public abstract void remove(VariantEvidenceSupport ves);
		public abstract Iterator<VariantEvidenceSupport> findAllOverlapping(BreakendSummary breakend);
		public abstract VariantEvidenceSupport findBestOverlapping(BreakendSummary breakend);
		protected void add(IntervalTree<List<VariantEvidenceSupport>> lookup, int start, int end, VariantEvidenceSupport ves) {
			Node<List<VariantEvidenceSupport>> node = lookup.find(start, end);
			if (node != null) {
				node.getValue().add(ves);
			} else {
				ArrayList<VariantEvidenceSupport> list = new ArrayList<>(2);
				list.add(ves);
				lookup.put(start, end, list);
			}
		}
		protected void remove(IntervalTree<List<VariantEvidenceSupport>> lookup, int start, int end, VariantEvidenceSupport ves) {
			Node<List<VariantEvidenceSupport>> node = lookup.find(start, end);
			if (node != null) {
				List<VariantEvidenceSupport> list = node.getValue();
				if (!list.remove(ves)) {
					String msg = String.format("Attempting to remove %s which does not exist on interval (%d, %d)", ves.location, start, end);
					throw new IllegalStateException(msg);
				}
				node.getValue().remove(ves);
				if (node.getValue().isEmpty()) {
					lookup.remove(start, end);
				}
			} else {
				String msg = String.format("Attempting to remove %s from non-existant interval (%d, %d)", ves.location, start, end);
				throw new IllegalStateException(msg);
			}
		}
		/**
		 * Gets the index of the IntervalTree lookup for this reference contig and direction  
		 */
		protected int getIndex(int referenceIndex, BreakendDirection dir) {
			return 2 * referenceIndex + (dir == BreakendDirection.Forward ? 0 : 1);
		}
		/**
		 * Creates an IntervalTree lookup for each reference contig and direction 
		 */
		protected List<IntervalTree<List<VariantEvidenceSupport>>> createByReferenceIndexDirectionLookup(int referenceSequenceCount) {
			return IntStream.range(0, referenceSequenceCount * 2)
					.mapToObj(i -> new IntervalTree<List<VariantEvidenceSupport>>())
					.collect(Collectors.toList());
		}
		protected VariantEvidenceSupport findBestOverlapping(BreakendSummary breakend, Iterator<VariantEvidenceSupport> it) {
			VariantEvidenceSupport best = null;
			while (it.hasNext()) {
				VariantEvidenceSupport v = it.next();
				assert(v.location.overlaps(breakend));
				if (best == null || ByScoreAscPositionDesc.compare(v, best) > 0) { 
					best = v;
				}
			}
			return best;
		}
	}
	/**
	 * Finds all variant calls overlapping the given evidence 
	 * @author Daniel Cameron
	 *
	 */
	private static class LocalOverlapLookup extends OverlapLookup {
		List<IntervalTree<List<VariantEvidenceSupport>>> localLookup;
		public LocalOverlapLookup(int referenceSequenceCount) {
			localLookup = createByReferenceIndexDirectionLookup(referenceSequenceCount);
		}
		public void add(VariantEvidenceSupport ves)
		{
			add(localLookup.get(getIndex(ves.location.referenceIndex, ves.location.direction)), ves.location.start, ves.location.end, ves);
		}
		public void remove(VariantEvidenceSupport ves)
		{
			remove(localLookup.get(getIndex(ves.location.referenceIndex, ves.location.direction)), ves.location.start, ves.location.end, ves);
		}
		public Iterator<VariantEvidenceSupport> findAllOverlapping(BreakendSummary breakend) {
			IntervalTree<List<VariantEvidenceSupport>> tree = localLookup.get(getIndex(breakend.referenceIndex, breakend.direction));
			Iterator<Node<List<VariantEvidenceSupport>>> nodeit = tree.overlappers(breakend.start, breakend.end);
			return new VariantEvidenceSupportNodeListIterator(nodeit, breakend);
		}
		public VariantEvidenceSupport findBestOverlapping(BreakendSummary breakend) {
			Iterator<VariantEvidenceSupport> it = findAllOverlapping(breakend);
			return findBestOverlapping(breakend, it);
		}
	}
	private static class VariantEvidenceSupportNodeListIterator extends AbstractIterator<VariantEvidenceSupport> {
		private final Iterator<Node<List<VariantEvidenceSupport>>> nodeIt;
		private final BreakendSummary breakend;
		private Iterator<VariantEvidenceSupport> listIt = Collections.emptyIterator();
		public VariantEvidenceSupportNodeListIterator(Iterator<Node<List<VariantEvidenceSupport>>> nodeIt, BreakendSummary breakend) {
			this.nodeIt = nodeIt;
			this.breakend = breakend;
		}
		@Override
		protected VariantEvidenceSupport computeNext() {
			while (nodeIt.hasNext() || listIt.hasNext()) {
				// advance through current list
				while (listIt.hasNext()) {
					VariantEvidenceSupport ves = listIt.next();
					if (ves.location.overlaps(breakend)) {
						return ves;
					}
				}
				if (nodeIt.hasNext()) {
					listIt = nodeIt.next().getValue().iterator();
				}
			}
			return endOfData();
		}
	}
	/**
	 * Finds the best overlapping variant call using a 1D interval tree on the remote breakpoint
	 * This should have better performance as, for repetitive sequence, the remote breakends
	 * are distributed across all repeats, but the local breakends all map to the same location
	 * (since we are doing a sequential traversal).
	 */
	private static class RemoteOverlapLookup extends OverlapLookup {
		List<IntervalTree<List<VariantEvidenceSupport>>> remoteLookup;
		List<RangeMap<Integer, VariantEvidenceSupport>> bestLocal;
		public RemoteOverlapLookup(int referenceSequenceCount) {
			this.remoteLookup = createByReferenceIndexDirectionLookup(referenceSequenceCount);
			this.bestLocal = IntStream.range(0, referenceSequenceCount * 2)
					.mapToObj(i -> TreeRangeMap.<Integer,VariantEvidenceSupport>create())
					.collect(Collectors.toList());
		}
		public void add(VariantEvidenceSupport ves) {
			assert(ves.location instanceof BreakpointSummary);
			BreakpointSummary location = (BreakpointSummary)ves.location;
			add(remoteLookup.get(getIndex(location.referenceIndex2, location.direction2)), location.start2, location.end2, ves);
			// need to add over the intervals in which we are the best
			RangeMap<Integer, VariantEvidenceSupport> rm = bestLocal.get(getIndex(location.referenceIndex, location.direction));
			RangeUtil.addWhereBest(rm, Range.closedOpen(location.start, location.end + 1), ves, ByScoreAscPositionDesc);
		}
		public void remove(VariantEvidenceSupport ves) {
			assert(ves.location instanceof BreakpointSummary);
			BreakpointSummary location = (BreakpointSummary)ves.location;
			remove(remoteLookup.get(getIndex(location.referenceIndex2, location.direction2)), location.start2, location.end2, ves);
			// we can remove all intervals before our end position as to be removed,
			// we need to have already added all the potential support for any variant
			// before our end position
			RangeMap<Integer, VariantEvidenceSupport> rm = bestLocal.get(getIndex(location.referenceIndex, location.direction));
			rm.remove(Range.lessThan(location.end + 1));
		}
		@Override
		public Iterator<VariantEvidenceSupport> findAllOverlapping(BreakendSummary breakend) {
			throw new NotImplementedException("RemoteOverlapLookup requires unique greedy evidence assignment");
		}
		public VariantEvidenceSupport findBestOverlapping(BreakpointSummary breakend) {
			IntervalTree<List<VariantEvidenceSupport>> tree = remoteLookup.get(getIndex(breakend.referenceIndex2, breakend.direction2));
			Iterator<Node<List<VariantEvidenceSupport>>> nodeit = tree.overlappers(breakend.start2, breakend.end2);
			Iterator<VariantEvidenceSupport> it = new VariantEvidenceSupportNodeListIterator(nodeit, breakend);
			return findBestOverlapping(breakend, it);
		}
		@Override
		public VariantEvidenceSupport findBestOverlapping(BreakendSummary breakend) {
			if (breakend instanceof BreakpointSummary) {
				return findBestOverlapping((BreakpointSummary)breakend);
			}
			// lookup bestLocal
			RangeMap<Integer, VariantEvidenceSupport> rm = bestLocal.get(getIndex(breakend.referenceIndex, breakend.direction));
			rm = rm.subRangeMap(Range.closedOpen(breakend.start, breakend.end + 1));
			float bestScore = -1;
			VariantEvidenceSupport best = null;
			for (VariantEvidenceSupport entry : rm.asMapOfRanges().values()) {
				if (entry.score > bestScore) {
					bestScore =  entry.score;
					best = entry;
				}
			}
			return best;
		}
	}
}
