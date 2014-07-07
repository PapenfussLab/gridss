package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import au.edu.wehi.idsv.AssemblyBuilder;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordUtil;
import au.edu.wehi.idsv.SoftClipEvidence;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;

/**
 * Debruijn graph implementation
 * @author Daniel Cameron
 *
 * @param <T>
 */
public abstract class DeBruijnGraphBase<T extends DeBruijnNodeBase> {
	public static final int MAX_QUAL_SCORE = 128 - 66;
	protected final Map<Long, T> kmers = Maps.newHashMap();
	protected final int k;
	protected final BreakendDirection direction;
	protected final ProcessingContext processContext;
	public DeBruijnGraphBase(ProcessingContext context, int k, BreakendDirection direction) {
		this.processContext = context;
		this.k = k;
		this.direction = direction;
	}
	public int getK() {
		return k;
	}
	public BreakendDirection getDirection() {
		return direction;
	}
	public void addEvidence(DirectedEvidence evidence) {
		if (evidence instanceof NonReferenceReadPair) {
			addEvidence((NonReferenceReadPair)evidence);
		} else if (evidence instanceof SoftClipEvidence) {
			addEvidence((SoftClipEvidence)evidence);
		} else {
			throw new RuntimeException(String.format("NYI: Unable to add %s evidence to de bruijn graph", evidence));
		}
	}
	protected void addEvidence(NonReferenceReadPair pair) {
		DeBruijnEvidence graphEvidence = DeBruijnEvidence.createRemoteReadEvidence(direction, k, pair);
		addEvidenceKmers(graphEvidence);
	}
	protected void addEvidence(SoftClipEvidence read) {
		DeBruijnEvidence graphEvidence = DeBruijnEvidence.createSoftClipEvidence(direction, k, read);
		addEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(NonReferenceReadPair pair) {
		DeBruijnEvidence graphEvidence = DeBruijnEvidence.createRemoteReadEvidence(direction, k, pair);
		removeEvidenceKmers(graphEvidence);
	}
	public void removeEvidence(SoftClipEvidence read) {
		DeBruijnEvidence graphEvidence = DeBruijnEvidence.createSoftClipEvidence(direction, k, read);
		removeEvidenceKmers(graphEvidence);
	}
	protected void addEvidenceKmers(DeBruijnEvidence evidence) {
		int readKmerOffset = 0;
		SAMRecord record = evidence.getSAMRecord();
		for (ReadKmer kmer : new ReadKmerIterable(k, record.getReadBases(), record.getBaseQualities(), evidence.isReversed(), evidence.isComplemented())) {
			if (evidence.isSkippedKmer(readKmerOffset)) {
				// do nothing with skipped kmers
			} else {
				addKmer(evidence, readKmerOffset, kmer);
			}
			readKmerOffset++;
		}
	}
	protected abstract T createEmptyNode();
	/**
	 * Adds the given kmer to the de bruijn graph
	 * @param evidence source evidence for this kmer
	 * @param readKmerOffset read offset of this kmer
	 * @param kmer kmer 
	 * @return de bruijn graph node for this kmer
	 */
	protected T addKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		T node = kmers.get(kmer.kmer);
		boolean newNode = node == null;
		if (newNode) {
			node = createEmptyNode();
			kmers.put(kmer.kmer, node);
		}
		node.add(evidence, readKmerOffset, kmer);
		if (newNode) {
			kmerAddedToGraph(evidence, readKmerOffset, kmer);
		}
		return node;
	}
	protected void removeEvidenceKmers(DeBruijnEvidence evidence) {
		int readKmerOffset = 0;
		SAMRecord record = evidence.getSAMRecord();
		for (ReadKmer kmer : new ReadKmerIterable(k, record.getReadBases(), record.getBaseQualities(), evidence.isReversed(), evidence.isComplemented())) {
			if (evidence.isSkippedKmer(readKmerOffset)) {
				// do nothing with skipped kmers
			} else {
				removeKmer(evidence, readKmerOffset, kmer);
			}
			readKmerOffset++;
		}
	}
	/**
	 * Removes the given kmer to the de bruijn graph
	 * @param evidence source evidence for this kmer
	 * @param readKmerOffset read offset of this kmer
	 * @param kmer kmer 
	 * @return de bruijn graph node for this kmer
	 */
	protected T removeKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		T node = kmers.get(kmer.kmer);
		if (node != null) {
			if (node.remove(evidence, readKmerOffset, kmer)) {
				kmers.remove(node);
				kmerRemovedFromGraph(evidence, readKmerOffset, kmer);
			}
		}
		return node;
	}
	protected void kmerAddedToGraph(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) { }
	protected void kmerRemovedFromGraph(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) { }
	/**
	 * Adjusts base qualities to be within valid FASTQ encoding range 
	 * @param bases base qualities to adjust
	 * @return 0-based phred-encodable base qualities
	 */
	protected byte[] rescaleBaseQualities(List<Integer> bases) {
		//Long largest = Collections.max(bases);
		//float scaleFactor = Math.min(1, MAX_QUAL_SCORE / (float)largest);
		byte[] result = new byte[bases.size()];
		for (int i = 0; i < result.length; i++) {
			//result[i] = (byte)(bases.get(i) * scaleFactor);
			result[i] = (byte)(bases.get(i) > MAX_QUAL_SCORE ? MAX_QUAL_SCORE : bases.get(i));
		}
		return result;
	}
	/**
	 * Gets the best kmer following the given kmer
	 * @param state kmer
	 * @param inclusionSet kmers must be in this set. Parameter is ignored if null
	 * @param exclusionSet kmers must not be in this set. Parameter is ignored if null
	 * @return next kmer, null if no valid kmer 
	 */
	protected Long greedyNextState(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		int best = Integer.MIN_VALUE;
		Long bestNode = null;
		for (Long next : nextStates(state, inclusionSet, exclusionSet)) {
			int weight = kmers.get(next).getWeight();
			if (weight > best) {
				bestNode = next;
				best = weight;
			}
		}
		return bestNode; 
	}
	protected List<Long> nextStates(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		List<Long> result = Lists.newArrayListWithCapacity(4);
		for (Long next : KmerEncodingHelper.nextStates(k, state)) {
			DeBruijnNodeBase node = kmers.get(next);
			if (node != null) {
				if ((inclusionSet == null || inclusionSet.contains(next)) && 
						(exclusionSet == null || !exclusionSet.contains(next))) {
					result.add(next);
				}
			}
		}
		return result; 
	}
	/**
	 * Gets the best kmer preceeding the given kmer
	 * @param state kmer
	 * @param inclusionSet kmers must be in this set. Parameter is ignored if null
	 * @param exclusionSet kmers must not be in this set. Parameter is ignored if null
	 * @return previous kmer, null if no valid kmer 
	 */
	protected Long greedyPrevState(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		int best = Integer.MIN_VALUE;
		Long bestNode = null;
		for (Long prev : prevStates(state, inclusionSet, exclusionSet)) {
			int weight = kmers.get(prev).getWeight();
			if (weight > best) {
				bestNode = prev;
				best = weight;
			}
		}
		return bestNode; 
	}
	protected List<Long> prevStates(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		List<Long> result = Lists.newArrayListWithCapacity(4);
		for (Long next : KmerEncodingHelper.prevStates(k, state)) {
			DeBruijnNodeBase node = kmers.get(next);
			if (node != null) {
				if ((inclusionSet == null || inclusionSet.contains(next)) && 
						(exclusionSet == null || !exclusionSet.contains(next))) {
					result.add(next);
				}
			}
		}
		return result; 
	}
	/**
	 * Base calls of contig
	 * @param path kmer contig
	 * @return base calls of a positive strand SAMRecord readout of contig
	 */
	protected byte[] getBaseCalls(List<Long> path) {
		int assemblyLength = path.size() + k - 1;
		byte[] bases = KmerEncodingHelper.encodedToPicardBases(k, path.get(0));
		bases = Arrays.copyOf(bases, assemblyLength);
		int offset = k - 1;
		for (Long node : path) {
			bases[offset] = KmerEncodingHelper.lastBaseEncodedToPicardBase(k, node);
			offset++;
		}
		if (direction == BreakendDirection.Backward) {
			ArrayUtils.reverse(bases);
		}
		return bases;
	}
	/**
	 * Base qualities of contig
	 * @param path kmer contig
	 * @return base qualities of a positive strand SAMRecord readout of contig
	 */
	protected byte[] getBaseQuals(List<Long> path) {
		List<Integer> qual = new ArrayList<Integer>(path.size());
		for (Long node : path) {
			// subtract # reads to adjust for the +1 qual introduced by ReadKmerIterable
			// to ensure positive node weights
			qual.add(this.kmers.get(node).getWeight() - this.kmers.get(node).getSupportingReads().size());
		}
		// pad out qualities to match the path length
		for (int i = 0; i < k - 1; i++) qual.add(qual.get(qual.size() - 1));
		byte[] quals = rescaleBaseQualities(qual);
		if (direction == BreakendDirection.Backward) {
			ArrayUtils.reverse(quals);
		}
		return quals;
	}
	protected Set<SAMRecord> getSupportingSAMRecords(Iterable<Long> path) {
		Set<SAMRecord> reads = Sets.newHashSet();
		for (Long kmer : path) {
			reads.addAll(kmers.get(kmer).getSupportingReads());
		}
		return reads;
	}
	/**
	 * Number of read bases supporting the given path
	 * @param path kmer contig
	 * @return number of read bases include in at least one kmer on the given kmer contig
	 */
	protected int getSAMRecordBaseCount(List<Long> path) {
		int readBaseCount = 0;
		Set<SAMRecord> lastNodeSupport = Sets.newHashSet();
		for (Long node : path) {
			for (SAMRecord read : this.kmers.get(node).getSupportingReads()) {
				if (lastNodeSupport.contains(read)) {
					readBaseCount++;
				} else {
					readBaseCount += k;
				}
			}
			lastNodeSupport = this.kmers.get(node).getSupportingReads();
		}
		return readBaseCount;
	}
	/**
	 * Gets the length of the longest read in the given contig
	 * @param supportingReads
	 * @return length of longest read
	 */
	protected int getMaxReadLength(List<Long> path) {
		int readLength = 0;
		for (Long kmer : path) {
			for (SAMRecord r : kmers.get(kmer).getSupportingReads()) {
				readLength = Math.max(readLength, r.getReadLength());
			}
		}
		return readLength;
	}
	/**
	 * Gets the length of the longest soft clip in the given contig
	 * @param supportingReads
	 * @return length of longest read
	 */
	protected int getMaxSoftClipLength(List<Long> path, int anchorLength) {
		int readLength = 0;
		int offset = 0;
		for (Long kmer : path) {
			// don't consider SC length of reference kmers - they may be SC support for a different breakend 
			if (offset >= anchorLength) {
				for (SAMRecord r : kmers.get(kmer).getSupportingReads()) {
					readLength = Math.max(readLength, direction == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipLength(r) : SAMRecordUtil.getStartSoftClipLength(r));
				}
			}
			offset++;
		}
		return readLength;
	}
	protected AssemblyBuilder debruijnContigAssembly(List<Long> path) {
		Set<SAMRecord> support = getSupportingSAMRecords(path);
		AssemblyBuilder builder = new AssemblyBuilder(processContext)
			.direction(direction)
			.assemblyBases(getBaseCalls(path))
			.assemblyBaseQuality(getBaseQuals(path))
			.assembledReadCount(support.size())
			.assembledBaseCount(getSAMRecordBaseCount(path))
			.longestSupportingRead(getMaxReadLength(path));
		return builder;
	}
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	protected Ordering<Long> ByKmerWeight = new Ordering<Long>() {
		public int compare(Long kmer1, Long kmer2) {
			return Ints.compare(getWeight(kmer1), getWeight(kmer2));
		}
		private int getWeight(long kmer) {
			int weight = 0;
			T node = kmers.get(kmer);
			if (node != null) weight = node.getWeight();
			return weight;
		}
	};
	@Override
	public String toString() {
		return toString(16);
	}
	public String toString(int maxNodesToPrint) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("De Bruijn graph: k=%d, %d kmers\n", k, kmers.size()));
		for (Long x : kmers.keySet()) {
			sb.append(printKmer(x));
			maxNodesToPrint--;
			if (maxNodesToPrint <= 0) break;
		}
		return sb.toString();
	}
	protected String printKmer(long x) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("%s(%d): %s",
				KmerEncodingHelper.toString(k, x),
				x,
				kmers.get(x)
				));
		sb.append(" from:{");
		for (Long y : KmerEncodingHelper.prevStates(k, x)) {
			DeBruijnNodeBase node = kmers.get(y);
			if (node != null) {
				sb.append(KmerEncodingHelper.toString(k, y));
				sb.append(',');
			}
		}
		sb.append("} to:{");
		for (Long y : KmerEncodingHelper.nextStates(k, x)) {
			DeBruijnNodeBase node = kmers.get(y);
			if (node != null) {
				sb.append(KmerEncodingHelper.toString(k, y));
				sb.append(',');
			}
		}
		sb.append("}\n");
		return sb.toString();
	}
	public String debugPrintPaths() {
		Map<Long, Integer> contigLookup = Maps.newHashMap();
		HashSet<Long> remaining = Sets.newHashSet(kmers.keySet());
		List<LinkedList<Long>> paths = Lists.newArrayList();
		int contig = 0;
		// enumerate the paths
		while (!remaining.isEmpty()) {
			contig++;
			LinkedList<Long> path = new LinkedList<Long>();
			path.add(remaining.iterator().next());
			remaining.remove(path.getFirst());
			contigLookup.put(path.getFirst(), contig);
			for(List<Long> adj = nextStates(path.getLast(), null, null); adj.size() == 1 && prevStates(adj.get(0), null, null).size() <= 1; adj = nextStates(path.getLast(), null, null)) {
				contigLookup.put(adj.get(0), contig);
				path.addLast(adj.get(0));
				remaining.remove(adj.get(0));
			}
			for(List<Long> adj = prevStates(path.getFirst(), null, null); adj.size() == 1 && nextStates(adj.get(0), null, null).size() <= 1; adj = prevStates(path.getFirst(), null, null)) {
				contigLookup.put(adj.get(0), contig);
				path.addFirst(adj.get(0));
				remaining.remove(adj.get(0));
			}
			paths.add(path);
		}
		return debugPrintPaths(paths, contigLookup);
	}
	protected String debugPrintPaths(List<LinkedList<Long>> paths, Map<Long, Integer> contigLookup) {
		StringBuilder sb = new StringBuilder();
		for (LinkedList<Long> path : paths) {
			sb.append(printPath(path, contigLookup));
		}
		return sb.toString();
	}
	protected String printPathAttributes(LinkedList<Long> path) {
		return "";
	}
	protected String printPath(LinkedList<Long> path, Map<Long, Integer> contigLookup) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("[%3d]\t", contigLookup.get(path.getFirst())));
		sb.append(new String(getBaseCalls(path)));
		int weight = 0;
		for (Long y : path) weight += kmers.get(y).getWeight();
		sb.append(String.format(" %d kmers %d weight", path.size(), weight));
		sb.append(" from:{");
		for (Long y : prevStates(path.getFirst(), null, null)) {
			sb.append(contigLookup.get(y));
			sb.append(',');
		}
		sb.append("} to:{");
		for (Long y : nextStates(path.getLast(), null, null)) {
			sb.append(contigLookup.get(y));
			sb.append(',');
		}
		sb.append("}");
		sb.append(printPathAttributes(path));
		sb.append('\n');
		return sb.toString();
	}
}
