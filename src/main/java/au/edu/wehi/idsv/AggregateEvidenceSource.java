package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.CloseableIterator;

public class AggregateEvidenceSource extends EvidenceSource implements Iterable<DirectedEvidence> {
	private List<SAMEvidenceSource> all;
	public AggregateEvidenceSource(ProcessingContext processContext, List<SAMEvidenceSource> reads, AssemblyEvidenceSource assemblies) {
		super(processContext, null);
		this.all = new ArrayList<>(reads);
		if (assemblies != null) {
			this.all.add(assemblies);
		}
	}
	@Override
	public CloseableIterator<DirectedEvidence> iterator() {
		return SAMEvidenceSource.mergedIterator(all, true);
	}
	public CloseableIterator<DirectedEvidence> iterator(QueryInterval intervals) {
		return SAMEvidenceSource.mergedIterator(all, intervals);
	}
	@Override
	public int getMaxConcordantFragmentSize() {
		return all.stream().mapToInt(source -> source.getMaxConcordantFragmentSize()).max().getAsInt();
	}
	@Override
	public int getMinConcordantFragmentSize() {
		return all.stream().mapToInt(source -> source.getMinConcordantFragmentSize()).min().getAsInt();
	}
	@Override
	public int getMaxReadLength() {
		return all.stream().mapToInt(source -> source.getMaxReadLength()).max().getAsInt();
	}
	@Override
	public int getMaxReadMappedLength() {
		return all.stream().mapToInt(source -> source.getMaxReadMappedLength()).max().getAsInt();
	}
}