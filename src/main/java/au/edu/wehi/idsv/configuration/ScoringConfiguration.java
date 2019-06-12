package au.edu.wehi.idsv.configuration;

import au.edu.wehi.idsv.model.*;
import org.apache.commons.configuration.Configuration;

import java.util.List;

public class ScoringConfiguration {
	public static final String CONFIGURATION_PREFIX = "scoring";
	public ScoringConfiguration(Configuration config) {
		config = config.subset(CONFIGURATION_PREFIX);
		switch (config.getString("model")) {
			case "EmpiricalLlr":
				model = new EmpiricalLlrModel();
				break;
			case "EmpiricalReferenceLikelihood":
				model = new EmpiricalReferenceLikelihoodModel();
				break;
			case "FastEmpiricalReferenceLikelihood":
				model = new FastEmpiricalReferenceLikelihoodModel();
				break;
			case "Mapq":
				model = new MapqModel();
				break;
			case "ReadCount":
				model = new ReadCountModel();
				break;
			default:
				throw new IllegalArgumentException(String.format("Unrecognised variant scoring model \"%s\"", config.getString("model")));
		}
		List<Object> exclude = config.getList("exclude");
		model = new ExclusionModel(model,
				exclude.contains("DiscordantPair"),
				exclude.contains("UnmappedMate"),
				exclude.contains("SplitRead"),
				exclude.contains("SoftClip"),
				exclude.contains("Indel"));
	}
	private VariantScoringModel model;
	public VariantScoringModel getModel() {
		return model;
	}
}
