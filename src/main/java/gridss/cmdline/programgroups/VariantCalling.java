package gridss.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class VariantCalling implements CommandLineProgramGroup {
	@Override
    public String getName() { return "VariantCalling"; }
    @Override
    public String getDescription() { return "Tools for calling variants."; }
}
