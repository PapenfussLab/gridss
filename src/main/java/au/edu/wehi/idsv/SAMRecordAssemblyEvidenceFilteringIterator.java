package au.edu.wehi.idsv;

import java.util.Iterator;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

public class SAMRecordAssemblyEvidenceFilteringIterator extends AbstractIterator<SAMRecordAssemblyEvidence> {
	private final ProcessingContext processContext;
	private final PeekingIterator<SAMRecordAssemblyEvidence> it;
	public SAMRecordAssemblyEvidenceFilteringIterator(
			ProcessingContext processContext,
			Iterator<SAMRecordAssemblyEvidence> it) {
		this.processContext = processContext;
		this.it = Iterators.peekingIterator(it); 
	}
	@Override
	protected SAMRecordAssemblyEvidence computeNext() {
		while (it.hasNext()) {
			SAMRecordAssemblyEvidence evidence = it.next();
			processContext.getAssemblyParameters().applyFilters(evidence);
			if (evidence.getBreakendSummary() instanceof BreakpointSummary) {
				for (VcfFilter breakpointFilter : processContext.getVariantCallingParameters().breakpointFilters((BreakpointSummary)evidence.getBreakendSummary())) {
					evidence.filterAssembly(breakpointFilter);
				}
			}
			if (!evidence.isAssemblyFiltered() || processContext.getAssemblyParameters().writeFilteredAssemblies) {
				return evidence;
			}
		}
		return endOfData();
	}
}
