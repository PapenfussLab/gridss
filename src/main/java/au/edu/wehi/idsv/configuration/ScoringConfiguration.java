package au.edu.wehi.idsv.configuration;

import java.util.List;

import org.apache.commons.configuration.Configuration;

import au.edu.wehi.idsv.model.EmpiricalLlrModel;
import au.edu.wehi.idsv.model.EmpiricalReferenceLikelihoodModel;
import au.edu.wehi.idsv.model.ExclusionModel;
import au.edu.wehi.idsv.model.FastEmpiricalReferenceLikelihoodModel;
import au.edu.wehi.idsv.model.MapqModel;
import au.edu.wehi.idsv.model.ReadCountModel;
import au.edu.wehi.idsv.model.VariantScoringModel;

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
