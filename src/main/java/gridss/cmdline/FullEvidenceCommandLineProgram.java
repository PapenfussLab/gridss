package gridss.cmdline;

import java.io.File;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
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
