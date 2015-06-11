package au.edu.wehi.idsv.util;

/**
 * Signals that an algorithm significantly exceeded the expected
 * computational time.
 *  
 * @author Daniel Cameron
 *
 */
public class AlgorithmRuntimeSafetyLimitExceededException extends RuntimeException {
	public AlgorithmRuntimeSafetyLimitExceededException() {
		super();
	}
	public AlgorithmRuntimeSafetyLimitExceededException(String message, Throwable ex) {
		super(message, ex);
	}
	public AlgorithmRuntimeSafetyLimitExceededException(String arg0) {
		super(arg0);
	}
	public AlgorithmRuntimeSafetyLimitExceededException(Throwable ex) {
		super(ex);
	}
	/**
	 * 
	 */
	private static final long serialVersionUID = 6139117344055780121L;

}
