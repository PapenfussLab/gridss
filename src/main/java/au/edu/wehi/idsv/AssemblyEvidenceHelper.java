package au.edu.wehi.idsv;

import java.util.Collection;

public class AssemblyEvidenceHelper {
	public static int getMapq(Collection<DirectedEvidence> support) {
		int mapq = 0;
		for (DirectedEvidence e : support) {
			mapq = Math.max(mapq, e.getLocalMapq());
		}
		return mapq;
	}
}
