package au.edu.wehi.idsv;

/**
 * Adds additional VCF attributes to the given breakend
 * @author cameron.d
 *
 */
public interface BreakendAnnotator {
	public abstract VariantContextDirectedEvidence annotate(VariantContextDirectedEvidence variant);

}