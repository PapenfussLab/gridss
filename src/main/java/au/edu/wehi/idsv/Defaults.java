package au.edu.wehi.idsv;



public class Defaults {
	public static final boolean SANITY_CHECK_ASSEMBLY_GRAPH;
	public static final boolean SANITY_CHECK_CLIQUE;
	public static final boolean SANITY_CHECK_ITERATORS;
	public static final boolean SANITY_CHECK_MEMOIZATION;
	public static final boolean SANITY_CHECK_MEMOIZATION_ALL_OPERATIONS;
	public static final boolean SINGLE_THREAD_LIBSSW;
	public static final boolean NO_LIBSSW;
	public static final boolean ASYNC_CACHE_REFERENCE;
	public static final boolean CACHE_REFERENCE;
	public static final boolean ATTEMPT_ASSEMBLY_RECOVERY;
	public static final boolean USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES;
	static {
		SANITY_CHECK_ASSEMBLY_GRAPH = Boolean.valueOf(System.getProperty("sanitycheck.assembly", "false"));
		SANITY_CHECK_CLIQUE = Boolean.valueOf(System.getProperty("sanitycheck.clique", "false"));
		SANITY_CHECK_ITERATORS = Boolean.valueOf(System.getProperty("sanitycheck.iterators", "false"));
		SANITY_CHECK_MEMOIZATION = Boolean.valueOf(System.getProperty("sanitycheck.memoization", "false"));
		SANITY_CHECK_MEMOIZATION_ALL_OPERATIONS = Boolean.valueOf(System.getProperty("sanitycheck.memoization.alloperations", "false"));
		SINGLE_THREAD_LIBSSW = Boolean.valueOf(System.getProperty("sswjni.sync", "false"));
		NO_LIBSSW = Boolean.valueOf(System.getProperty("sswjni.disable", "false"));
		ASYNC_CACHE_REFERENCE = !Boolean.valueOf(System.getProperty("reference.loading.sync", "false"));
		CACHE_REFERENCE = !Boolean.valueOf(System.getProperty("reference.cache", "true"));
		ATTEMPT_ASSEMBLY_RECOVERY = Boolean.valueOf(System.getProperty("assembly.recover", "true"));
		USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES = Boolean.valueOf(System.getProperty("assembly.optimised_data_structures", "true"));
	}
}