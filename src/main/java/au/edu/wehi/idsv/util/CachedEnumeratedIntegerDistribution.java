package au.edu.wehi.idsv.util;

import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.hash.TIntDoubleHashMap;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.NotANumberException;
import org.apache.commons.math3.exception.NotFiniteNumberException;
import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Caches results so the underlying distribution is not converted
 * to/from an array for every function call
 * @author cameron.d
 *
 */
public class CachedEnumeratedIntegerDistribution extends EnumeratedIntegerDistribution {
    /**
	 * 
	 */
	private static final long serialVersionUID = 126014755178519095L;
	public CachedEnumeratedIntegerDistribution(final int[] singletons, final double[] probabilities)
    throws DimensionMismatchException, NotPositiveException, MathArithmeticException,
           NotFiniteNumberException, NotANumberException{
        super(singletons, probabilities);
    }
    public CachedEnumeratedIntegerDistribution(final RandomGenerator rng, final int[] singletons, final double[] probabilities)
        throws DimensionMismatchException, NotPositiveException, MathArithmeticException,
                NotFiniteNumberException, NotANumberException {
    	super(rng, singletons, probabilities);
    }
    private TIntDoubleMap cacheProbability = new TIntDoubleHashMap();
    public synchronized double probability(final int x) {
    	if (!cacheProbability.containsKey(x) && x >= getSupportLowerBound() && x <= getSupportUpperBound()) {
    		cacheProbability.put(x, super.probability(x));
    	}
		return cacheProbability.get(x);
    }
    private TIntDoubleMap cacheCumulativeProbability = new TIntDoubleHashMap();
    public synchronized double cumulativeProbability(final int x) {
    	if (!cacheCumulativeProbability.containsKey(x) && x >= getSupportLowerBound() && x <= getSupportUpperBound()) {
    		cacheCumulativeProbability.put(x, super.cumulativeProbability(x));
    	}
		return cacheCumulativeProbability.get(x);
    }
    private Double cacheNumericalMean;
    public double getNumericalMean() {
    	if (cacheNumericalMean == null) {
    		cacheNumericalMean = super.getNumericalMean();
    	}
    	return cacheNumericalMean;
    }
    private Double cacheNumericalVariance;
    public double getNumericalVariance() {
    	if (cacheNumericalVariance == null) {
    		cacheNumericalVariance = super.getNumericalVariance();
    	}
    	return cacheNumericalVariance;
    }
    private Integer cacheSupportLowerBound;
    public int getSupportLowerBound() {
    	if (cacheSupportLowerBound == null) {
    		cacheSupportLowerBound = super.getSupportLowerBound();
    	}
    	return cacheSupportLowerBound;
    }
    private Integer cacheSupportUpperBound;
    public int getSupportUpperBound() {
    	if (cacheSupportUpperBound == null) {
    		cacheSupportUpperBound = super.getSupportUpperBound();
    	}
    	return cacheSupportUpperBound;
    }
}
