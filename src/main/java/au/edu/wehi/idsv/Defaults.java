package au.edu.wehi.idsv;



public class Defaults {
	public static final boolean SANITY_CHECK_DE_BRUIJN;
	public static final boolean SANITY_CHECK_CLIQUE;
	public static final boolean SANITY_CHECK_ITERATORS;
	public static final boolean SINGLE_THREAD_LIBSSW;
	static {
		SANITY_CHECK_DE_BRUIJN = Boolean.valueOf(System.getProperty("sanitycheck.debruijn", "false"));
		SANITY_CHECK_CLIQUE = Boolean.valueOf(System.getProperty("sanitycheck.clique", "false"));
		SANITY_CHECK_ITERATORS = Boolean.valueOf(System.getProperty("sanitycheck.iterators", "false"));
		SINGLE_THREAD_LIBSSW = Boolean.valueOf(System.getProperty("sswjni.sync", "false"));
	}
}