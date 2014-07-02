package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.SAMRecord;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Debruijn graph implementation
 * @author Daniel Cameron
 *
 * @param <T>
 */
public abstract class DeBruijnGraphBase<T extends DeBruijnNodeBase> {
	public static final int MAX_QUAL_SCORE = 128 - 66;
	private final Map<Long, T> kmers = Maps.newHashMap();
	private final int k;
	private final BreakendDirection direction;
	public DeBruijnGraphBase(int k, BreakendDirection direction) {
		this.k = k;
		this.direction = direction;
	}
	public void addEvidence(NonReferenceReadPair pair) {
		DeBruijnEvidence graphEvidence = DeBruijnEvidence.createRemoteReadEvidence(direction, k, pair);
		addEvidenceKmers(graphEvidence);
	}
	public void addEvidence(SoftClipEvidence read) {
		DeBruijnEvidence graphEvidence = DeBruijnEvidence.createSoftClipEvidence(direction, k, read);
		addEvidenceKmers(graphEvidence);
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
	protected void addKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		
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
	protected void removeKmer(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		T node = kmers.get(kmer.kmer);
		if (node != null) {
			if (node.remove(evidence, readKmerOffset, kmer)) {
				kmers.remove(node);
				kmerRemovedFromGraph(evidence, readKmerOffset, kmer);
			}
		}
	}
	protected void kmerAddedToGraph(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) { }
	protected void kmerRemovedFromGraph(DeBruijnEvidence evidence, int readKmerOffset, ReadKmer kmer) { }
	protected Set<SAMRecord> getPathSupportingReads(LinkedList<Long> path) {
		Set<SAMRecord> reads = Sets.newHashSet();
		for (Long kmer : path) {
			reads.addAll(kmers.get(kmer).getSupportingReads());
		}
		return reads;
	}
	/**
	 * Adjusts base qualities to be within valid FASTQ encoding range 
	 * @param bases base qualities to adjust
	 * @return 0-based phred-encodable base qualities
	 */
	protected byte[] rescaleBaseQualities(List<Long> bases) {
		//Long largest = Collections.max(bases);
		//float scaleFactor = Math.min(1, MAX_QUAL_SCORE / (float)largest);
		byte[] result = new byte[bases.size()];
		for (int i = 0; i < result.length; i++) {
			//result[i] = (byte)(bases.get(i) * scaleFactor);
			result[i] = (byte)(bases.get(i) > MAX_QUAL_SCORE ? MAX_QUAL_SCORE : bases.get(i));
		}
		return result;
	}
	protected LinkedList<Long> greedyTraverseForward(Long start) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = new HashSet<Long>();
		path.add(start);
		visited.add(start);
		for (Long node = greedyNextState(start, visited); node != null; node = greedyNextState(node, visited)) {
			path.addLast(node);
			visited.add(node);
		}
		return path;
	}
	protected LinkedList<Long> greedyTraverseBackward(Long start) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = new HashSet<Long>();
		path.add(start);
		visited.add(start);
		for (Long node = greedyPrevState(start, visited); node != null; node = greedyPrevState(node, visited)) {
			path.addFirst(node);
			visited.add(node);
		}
		return path;
	}
	protected Long greedyNextState(Long state, Set<Long> visited) {
		long best = -1;
		Long bestNode = null;
		for (Long next : KmerEncodingHelper.nextStates(state, k)) {
			DeBruijnNodeBase node = kmers.get(next);
			if (node != null && node.getWeight() > best && !visited.contains(next)) {
				bestNode = next;
				best = node.getWeight();
			}
		}
		return bestNode; 
	}
	protected Long greedyPrevState(Long state, Set<Long> visited) {
		long best = -1;
		Long bestNode = null;
		for (Long next : KmerEncodingHelper.prevStates(state, k)) {
			DeBruijnNodeBase node = kmers.get(next);
			if (node != null && node.getWeight() > best && !visited.contains(next)) {
				bestNode = next;
				best = node.getWeight();
			}
		}
		return bestNode; 
	}
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("De Bruijn graph: k=%d, %d kmers\n", k, kmers.size()));
		int max = 10;
		for (Long x : kmers.keySet()) {
			sb.append(String.format("%s(%d): %d weight from %d reads",
					KmerEncodingHelper.toString(k, x),
					x,
					kmers.get(x).getWeight(),
					kmers.get(x).getSupportingReads().size()
					));
			sb.append(" from:{");
			for (Long y : KmerEncodingHelper.prevStates(x, k)) {
				DeBruijnNodeBase node = kmers.get(y);
				if (node != null) {
					sb.append(KmerEncodingHelper.toString(k, y));
					sb.append(',');
				}
			}
			sb.append("} to:{");
			for (Long y : KmerEncodingHelper.nextStates(x, k)) {
				DeBruijnNodeBase node = kmers.get(y);
				if (node != null) {
					sb.append(KmerEncodingHelper.toString(k, y));
					sb.append(',');
				}
			}
			sb.append("}\n");
			max--;
			if (max <= 0) break;
		}
		return sb.toString();
	}
}
