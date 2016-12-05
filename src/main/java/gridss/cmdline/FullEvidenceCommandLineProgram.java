package gridss.cmdline;

import java.io.File;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.validation.OrderAssertingIterator;
import au.edu.wehi.idsv.validation.PairedEvidenceTracker;
import htsjdk.samtools.util.CloseableIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;

public abstract class FullEvidenceCommandLineProgram extends MultipleSamFileCommandLineProgram {
	@Option(doc="Breakend assemblies which have undergone split read identification", optional=false)
	public File ASSEMBLY;
	private final boolean requireAssembly;
	public FullEvidenceCommandLineProgram() {
		this(true);
	}
	public FullEvidenceCommandLineProgram(boolean requireAssembly) {
		this.requireAssembly = requireAssembly;
	}
	private AssemblyEvidenceSource assemblyEvidenceSource;
	public AssemblyEvidenceSource getAssemblySource() {
		if (assemblyEvidenceSource == null) {
			assemblyEvidenceSource = new AssemblyEvidenceSource(getContext(), getSamEvidenceSources(), ASSEMBLY); 
		}
		return assemblyEvidenceSource;
	}
	public void setAssemblySource(AssemblyEvidenceSource source) {
		assemblyEvidenceSource = source;
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
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = assemblyCustomCommandLineValidation();
		if (val != null) return val;
		return super.customCommandLineValidation();
	}
	public String[] assemblyCustomCommandLineValidation() {
		if (requireAssembly && ASSEMBLY != null && !ASSEMBLY.exists()) {
			return new String[] { "Missing ASSEMBLY file" };
		}
    	return null;
	}
	@Override
	public void copyInputs(CommandLineProgram cmd) {
		super.copyInputs(cmd);
		if (cmd instanceof FullEvidenceCommandLineProgram) {
			FullEvidenceCommandLineProgram prog = (FullEvidenceCommandLineProgram) cmd;
			prog.ASSEMBLY = ASSEMBLY;
			prog.assemblyEvidenceSource = assemblyEvidenceSource;
		}
	}
}
