package gridss.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class Metrics implements CommandLineProgramGroup {
	@Override
    public String getName() { return "Data Cleaning"; }
    @Override
    public String getDescription() { return "Tools for fixing data errors."; }
}
