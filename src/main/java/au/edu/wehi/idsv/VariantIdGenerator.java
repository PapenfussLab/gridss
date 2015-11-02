package au.edu.wehi.idsv;

/**
 * Generates unique variant identifiers
 * 
 * @author Daniel Cameron
 *
 */
public interface VariantIdGenerator {
	/**
	 * Generates an identifier for the given structural variant
	 * 
	 * @param breakpoint breakpoint with low coordinates local
	 * @return unique variant identifier
	 */
	String generate(BreakpointSummary breakpoint);
}
