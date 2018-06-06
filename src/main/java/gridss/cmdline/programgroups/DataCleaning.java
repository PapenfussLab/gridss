package gridss.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class DataCleaning implements CommandLineProgramGroup {
	@Override
    public String getName() { return "Metrics"; }
    @Override
    public String getDescription() { return "Tools for generating summary statistics."; }
}
