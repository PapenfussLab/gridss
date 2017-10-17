package gridss.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class Benchmarking implements CommandLineProgramGroup {
	@Override
    public String getName() { return "Benchmarking"; }
    @Override
    public String getDescription() { return "Tools for creating simulated variants and benchmarking variant callers."; }
}
