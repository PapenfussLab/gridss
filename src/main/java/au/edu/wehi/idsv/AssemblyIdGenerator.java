package au.edu.wehi.idsv;

/**
 * Generates unique variant identifiers
 * 
 * @author Daniel Cameron
 *
 */
public interface AssemblyIdGenerator {
	/**
	 * Generates an identifier for the given structural variant assembly
	 * 
	 */
	String generate(BreakendSummary breakpoint, byte[] baseCalls, int startAnchoredBaseCount, int endAnchoredBaseCount);
}
