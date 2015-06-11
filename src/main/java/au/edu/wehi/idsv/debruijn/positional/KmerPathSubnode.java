package au.edu.wehi.idsv.debruijn.positional;

import java.util.ArrayList;
import java.util.List;

public class KmerPathSubnode {
	private final KmerPathNode n;
	private final int start;
	private final int end;
	public KmerPathSubnode(KmerPathNode n, int start, int end) {
		assert(n.startPosition(0) <= start);
		assert(n.endPosition(0) >= end);
		assert(end - start >= 0);
		this.n = n;
		this.start = start;
		this.end = end;
	}
	public KmerPathNode node() { return n; }
	public int firstKmerStartPosition() { return start; }
	public int firstKmerEndPosition() { return end; }
	public List<KmerPathSubnode> next() {
		List<KmerPathSubnode> adj = new ArrayList<KmerPathSubnode>(n.next().size());
		int targetStart = start + n.length() + 1;
		int targetEnd = end + n.length() + 1;
		for (KmerPathNode pn : n.next()) {
			int pnStart = pn.startPosition(0);
			int pnEnd = pn.endPosition(0);
			// since next() is sorted, we only need to process the neighbours overlapping our interval
			if (pnEnd < targetStart) {
				continue;
			} else if (pnStart > targetEnd) {
				break;
			} else {
				adj.add(new KmerPathSubnode(pn,
						Math.max(targetStart, pnStart),
						Math.min(targetEnd, pnEnd)));
			}
		}
		return adj;
	}
	public List<KmerPathSubnode> prev() {
		List<KmerPathSubnode> adj = new ArrayList<KmerPathSubnode>(n.prev().size());
		int targetStart = start - 1;
		int targetEnd = end - 1;
		for (KmerPathNode pn : n.prev()) {
			int pnStart = pn.startPosition();
			int pnEnd = pn.endPosition();
			if (pnEnd < targetStart) {
				continue;
			} else if (pnStart > targetEnd) {
				break;
			} else {
				adj.add(new KmerPathSubnode(pn,
						Math.max(targetStart, pnStart) - pn.length(),
						Math.min(targetEnd, pnEnd) - pn.length()));
			}
		}
		return adj;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + ((n == null) ? 0 : n.hashCode());
		result = prime * result + start;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		KmerPathSubnode other = (KmerPathSubnode) obj;
		if (end != other.end)
			return false;
		if (n == null) {
			if (other.n != null)
				return false;
		} else if (!n.equals(other.n))
			return false;
		if (start != other.start)
			return false;
		return true;
	}
}
