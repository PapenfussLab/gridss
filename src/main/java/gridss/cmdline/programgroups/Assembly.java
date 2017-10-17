package gridss.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class Assembly implements CommandLineProgramGroup {
	@Override
    public String getName() { return "Assembly"; }
    @Override
    public String getDescription() { return "Tools for performing assembly."; }
}
