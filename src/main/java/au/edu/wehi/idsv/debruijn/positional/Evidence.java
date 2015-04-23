package au.edu.wehi.idsv.debruijn.positional;

import java.util.List;

public abstract class Evidence {
	List<KmerSupportNode> support;
	String id;
	abstract int start();
	// first position is: [start - error, start + error]
	abstract int errorWidth(); 
	abstract boolean isAnchored(int offset);
}
