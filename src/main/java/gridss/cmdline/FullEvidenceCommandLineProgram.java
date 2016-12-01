package gridss.cmdline;

import java.io.File;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.EvidenceSource;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import htsjdk.samtools.util.CloseableIterator;
import picard.cmdline.Option;

public abstract class FullEvidenceCommandLineProgram extends MultipleSamFileCommandLineProgram {
	@Option(doc="Breakend assemblies which have undergone split read identification")
	public File ASSEMBLY;
	private AssemblyEvidenceSource assemblyEvidenceSource;
	public AssemblyEvidenceSource getAssemblySource() {
		if (assemblyEvidenceSource == null) {
			assemblyEvidenceSource = new AssemblyEvidenceSource(getContext(), getSamEvidenceSources(), ASSEMBLY); 
		}
		return assemblyEvidenceSource;
	}
	public CloseableIterator<DirectedEvidence> getEvidenceIterator() {
		CloseableIterator<DirectedEvidence> evidenceIt;
		boolean assemblyOnly = getContext().getVariantCallingParameters().callOnlyAssemblies;
		if (assemblyOnly) {
			evidenceIt = SAMEvidenceSource.mergedIterator(ImmutableList.of(getAssemblySource()));
		} else {
			evidenceIt = SAMEvidenceSource.mergedIterator(ImmutableList.<SAMEvidenceSource>builder().addAll(getSamEvidenceSources()).add(getAssemblySource()).build());
		}
		if (Defaults.SANITY_CHECK_ITERATORS) {
			evidenceIt = new AutoClosingIterator<>(
					new PairedEvidenceTracker<>("Evidence",
							new OrderAssertingIterator<>(evidenceIt, DirectedEvidenceOrder.ByNatural)), evidenceIt);
		}
		return evidenceIt;
	}
	/**
	 * Maximum distance between the SAM alignment location of evidence, and the extrema of the
	 * breakend position supported by that evidence. 
	 * @return maximum order-of-order distance between evidence ordered by SAM alignment position and the breakend start position 
	 */
	public int maximumWindowSize() {
		int maxSize = 0;
		for (EvidenceSource source : getSamEvidenceSources()) {
			SAMEvidenceSource samSource = (SAMEvidenceSource)source;
			maxSize = Math.max(samSource.getMaxConcordantFragmentSize(), Math.max(samSource.getMaxReadLength(), samSource.getMaxReadMappedLength()));
		}
		maxSize = Math.max(maxSize, getAssemblySource().getMaxAssemblyLength()) + getContext().getVariantCallingParameters().maxBreakendHomologyLength;
		return maxSize + 2 * (getContext().getVariantCallingParameters().breakendMargin + 1);
	}
	@Override
	protected String[] customCommandLineValidation() {
		if (ASSEMBLY == null || !ASSEMBLY.exists()) {
            return new String[]{"Missing ASSEMBLY file"};
        }
		return super.customCommandLineValidation();
	}
}
