package gridss.cmdline;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import com.google.common.collect.ImmutableList;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public abstract class FullEvidenceCommandLineProgram extends MultipleSamFileCommandLineProgram {
	@Argument(doc="Breakend assemblies which have undergone split read identification", optional=false)
	public List<File> ASSEMBLY;
	private final boolean requireAssembly;
	public FullEvidenceCommandLineProgram() {
		this(true);
	}
	public FullEvidenceCommandLineProgram(boolean requireAssembly) {
		this.requireAssembly = requireAssembly;
	}
	private List<AssemblyEvidenceSource> assemblyEvidenceSource;
	public List<AssemblyEvidenceSource> getAssemblySource() {
		if (assemblyEvidenceSource == null) {
			assemblyEvidenceSource = new ArrayList<>();
			for (File ass : ASSEMBLY) {
				assemblyEvidenceSource.add(new AssemblyEvidenceSource(getContext(), getSamEvidenceSources(), ass));
			}
			AssemblyEvidenceSource.validateAllCategoriesAssembled(getContext(), assemblyEvidenceSource);
		}
		return assemblyEvidenceSource;
	}
	public void setAssemblySource(List<AssemblyEvidenceSource> source) {
		assemblyEvidenceSource = source;
    }
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = assemblyCustomCommandLineValidation();
		if (val != null) return val;
		return super.customCommandLineValidation();
	}
	public String[] assemblyCustomCommandLineValidation() {
		if (requireAssembly) {
			if (ASSEMBLY == null || ASSEMBLY.size() == 0) {
				return new String[]{"Missing ASSEMBLY file"};
			}
			for (int i = 0; i < ASSEMBLY.size(); i++) {
				File left = ASSEMBLY.get(i);
				for (int j = i + 1; j < ASSEMBLY.size(); j++) {
					File right = ASSEMBLY.get(j);
					if (left.getAbsolutePath().equals(right.getAbsolutePath())) {
						return new String[] {"Multiple ASSEMBLY files with same name."};
					}
					if (WORKING_DIR != null && left.getName().equals(right.getName())) {
						return new String[] {"ASSEMBLY files must have unique names."};
					}
				}
			}
			String[] err = getWorkingDirectoryFilenameCollisions(ImmutableList.of(ASSEMBLY, INPUT), WORKING_DIR);
			if (err != null) {
				return err;
			}
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
