package au.edu.wehi.idsv;

/**
 * Exception indicating that something that should never happen, happened.
 * @author Daniel Cameron
 *
 */
public class SanityCheckFailureException extends RuntimeException {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public SanityCheckFailureException(String msg) {
		super(msg);
	}

}
