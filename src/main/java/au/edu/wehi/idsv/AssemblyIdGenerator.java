package au.edu.wehi.idsv;

/**
 * Generates unique variant identifiers
 * 
 * @author cameron.d
 *
 */
public interface AssemblyIdGenerator {
	/**
	 * Generates an identifier for the given structural variant assembly
	 * 
	 */
	String generate(BreakendSummary breakpoint, byte[] baseCalls, int anchoredBaseCount);
}
