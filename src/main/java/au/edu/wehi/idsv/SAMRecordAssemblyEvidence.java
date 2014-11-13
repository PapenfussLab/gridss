package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.Map;
import java.util.Set;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;

public class SAMRecordAssemblyEvidence implements AssemblyEvidence {
	private SAMRecord record;
	private AssemblyEvidenceSource source;
	private BreakendSummary breakend;
	public SAMRecordAssemblyEvidence(BreakendSummary breakend,
			AssemblyEvidenceSource source,
			Set<DirectedEvidence> evidence,
			int anchoredBases,
			byte[] baseCalls,
			byte[] baseQuals,
			Map<VcfAttributes, Integer> attributes
			) {
		this.breakend = breakend;
		this.source = source;
		this.record = new SAMRecord();
	}
}
