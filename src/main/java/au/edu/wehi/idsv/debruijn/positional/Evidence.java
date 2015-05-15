package au.edu.wehi.idsv.debruijn.positional;

import java.util.BitSet;

import au.edu.wehi.idsv.debruijn.PackedReadKmerList;


public class Evidence extends PackedReadKmerList {
	private String id;
	private BitSet anchor;
	private int start;
	private int error;
	public KmerSupportNode node(int offset) {
		return new KmerSupportNode(this, offset);
	}
	public String evidenceId() { return id; }
	public int startPosition() { 
		return start;
	}
	public int errorWidth() {
		return error;
	}
	public boolean isAnchored(int offset) {
		return anchor.get(offset);
	}
	public Evidence(
			String evidenceId,
			int start,
			int error,
			int k, byte[] bases, byte[] qual, boolean reverse, boolean complement) {
		super(k, bases, qual, reverse, complement);
		assert(evidenceId != null);
		this.id = evidenceId;
		this.start = start;
		this.error = error;
	}
	@Override
	public int hashCode() {
		return id.hashCode();
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Evidence other = (Evidence) obj;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}
}
