package au.edu.wehi.idsv;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import htsjdk.samtools.util.Log;

/**
 * Performs GC bias adjustments according to the given file.
 * 
 * This file must contain two tab-seperated columns without any header lines.
 * The first column must contain the GC percentage for adjustment, the second the adjustment multipler.
 * The first column should contain the integer values 0-100 inclusive, the second a floating point adjustment multiplier.
 * @author Daniel Cameron
 *
 */
public class PrecomputedGcBiasAdjuster implements GcBiasAdjuster {
	private static final Log log = Log.getInstance(PrecomputedGcBiasAdjuster.class);
	private final double[] adjustment;
	public PrecomputedGcBiasAdjuster(File tsv) throws IOException {
		this.adjustment = load(tsv);
	}
	public static double[] load(File tsv) throws IOException {
		double[] adjustment = new double[101];
		Arrays.fill(adjustment, -1); // sentinal
		try (BufferedReader br = new BufferedReader(new FileReader(tsv))) {
		    String line = br.readLine();
		    while (line != null) {
		    	String[] split = line.split("\\s+");
		    	int gcPercentage = Integer.parseInt(split[0]);
		    	double value = Double.parseDouble(split[1]);
		    	if (gcPercentage <  0 || gcPercentage > 100) {
		    		String msg = String.format("Invalid gc percentage %s. Ignoring.", split[0]);
		    		log.error(msg);
		    		throw new IllegalArgumentException(msg);
		    	} else {
		    		if (value < 0) {
		    			String msg = String.format("GC percentage adjustment multipler for %d is %s. Negative multipliers are invalid.", gcPercentage, split[1]);
		    			log.error(msg);
		    			throw new IllegalArgumentException(msg);
		    		} else {
		    			adjustment[gcPercentage] = value;
		    		}
		    	}
    			line = br.readLine();
		    }
		}
		for (int i = 0; i <= 100; i++) {
			if (adjustment[i] == -1) {
				String msg = String.format("GC percentage adjustment multipler missing for gc percentage %d", i);
				log.error(msg);
    			throw new IllegalArgumentException(msg);
			}
		}
		return adjustment;
	}
	@Override
	public double adjustmentMultiplier(int gcPercentage) {
		if (gcPercentage < 0) throw new IllegalArgumentException("GC percentage cannot be negative");
		if (gcPercentage > 100) throw new IllegalArgumentException("GC percentage cannot be greater than 100");
		return adjustment[gcPercentage];
	}

}
