package au.edu.wehi.idsv.debruijn.positional;



public class KmerSupportNode extends KmerNode {
	public long getKmer() { return kmer; }
	public int getStartPosition() { return evidence.start() - evidence.errorWidth() + offset; }
	public int getEndPosition() { return evidence.start() + evidence.errorWidth() + offset; }
	public int getWeight() { return weight; }
	public boolean isReference() { return evidence.isAnchored(offset); }
	private final long kmer;
	private final int weight;
	private final int offset;
	private final Evidence evidence;
	public KmerSupportNode(Evidence evidence, int offset, long kmer, int weight) {
		this.evidence = evidence;
		this.offset = offset;
		this.kmer = kmer;
		this.weight = weight;
	}
}
/*
	public KmerNode(int offset, Evidence e) {
		this.kmer = e.kmers[offset];
		this.startPosition = e.start() - e.errorWidth() + offset;
		this.endPosition = e.start() + e.errorWidth() + offset;
		this.weight = e.weight[offset];
		this.support = new ArrayList<KmerNodeSupport>(4);
		this.support.add(new KmerNodeSupport(offset, e));
	}
	public KmerNode(long kmer, int startPosition, int endPosition, int weight, Collection<KmerNodeSupport> support) {
		this.kmer = kmer;
		this.startPosition = startPosition;
		this.endPosition = endPosition;
		this.weight = weight;
		this.support = new ArrayList<KmerNodeSupport>(support.size() + 4);
		this.support.addAll(support);
	}
	public KmerNode splitStartingAt(int startPosition) {
		assert(this.startPosition < startPosition);
		assert(this.endPosition >= startPosition);
		this.endPosition = startPosition;
		KmerNode following = new KmerNode(
			this.kmer,
			startPosition,
			this.endPosition,
			this.weight,
			this.support);
		return following;
	}
	public KmerNode subInterval(int startPosition, int endPosition) {
		assert(this.startPosition <= startPosition);
		assert(this.endPosition >= endPosition);
		if (this.startPosition == startPosition && this.endPosition == endPosition) {
			return this;
		} else {
			return new KmerNode(
				this.kmer,
				startPosition,
				endPosition,
				this.weight,
				this.support);
		}
	}
	public void add(int weight, Collection<KmerNodeSupport> support) {
		this.weight += weight;
		this.support.addAll(support);
	}
	public void add(KmerNode node) {
		assert(node.kmer == kmer);
		assert(node.startPosition == startPosition);
		assert(node.endPosition == endPosition);
		add(node.weight, node.support);
	}
	public boolean isAnchored() {
		for (KmerNodeSupport s : support) {
			if (s.evidence.isAnchored(s.offset)) {
				return true;
			}
		}
		return false;
	}
	public boolean removeSupport(Set<Evidence> toRemove) {
		// TODO: better than O(n^2) removal
		for (int i = support.size(); i >= 0 ; i++) {
			if (toRemove.contains(support.get(i))) {
				weight -= support.get(i).evidence.weight[support.get(i).offset];
				support.remove(i);
			}
		}
		return support.isEmpty();
	}
	public Collection<KmerNodeSupport> getSupport() { return support; }
	//PathNode containedIn;
	//Node equivalent; // path collapsed
	// baseDiff(equivalent.kmer, kmer) < threshold (threshold applies over entire path)
	// equivalent.position = node.position
	public boolean assertInvariants() {
		assert(support != null);
		assert(support.size() > 0);
		assert(weight > 0);
		assert(startPosition >= 0);
		assert(endPosition >= startPosition);
		int supportWeight = 0;
		for (KmerNodeSupport s : support) {
			assert(s.evidence != null);
			assert(s.offset < s.evidence.kmers.length);
			assert(s.offset >= 0);
			assert(s.evidence.kmers[s.offset] == kmer);
			supportWeight += s.evidence.weight[s.offset];
		}
		assert(supportWeight == weight);
		return true;
	}
}
*/