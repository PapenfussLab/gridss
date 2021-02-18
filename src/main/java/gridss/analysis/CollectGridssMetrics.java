package gridss.analysis;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.NotImplementedException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.analysis.CollectMultipleMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.argumentcollections.RequiredOutputArgumentCollection;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Class that is designed to instantiate and execute multiple metrics programs that extend
 * SinglePassSamProgram while making only a single pass through the SAM file and supplying
 * each program with the records as it goes.
 *
 */
@CommandLineProgramProperties(
        summary = "Extension of picard.CollectMultipleMetrics to include GRIDSS metrics. "
        		+ "Additional metrics are CollectCigarMetrics, CollectIdsvMetrics, CollectTagMetrics, and CollectMapqMetrics.",
        oneLineSummary = "A \"meta-metrics\" calculating program that produces multiple metrics for the provided SAM/BAM",
        programGroup = gridss.cmdline.programgroups.Metrics.class
)
public class CollectGridssMetrics extends CollectMultipleMetrics {
    public static enum GridssProgram {
    	CollectCigarMetrics,
    	CollectMapqMetrics,
    	CollectTagMetrics,
        CollectIdsvMetrics,
        ReportThresholdCoverage,
    }
    @Argument(doc = "Set of gridss metrics programs to apply during the pass through the SAM file.")
    public Set<GridssProgram> GRIDSS_PROGRAM = new LinkedHashSet<>(Arrays.asList(GridssProgram.values()));
    
    @Argument(doc = "Threshold coverage to report for ReportThresholdCoverage.", optional=true)
    public Integer THRESHOLD_COVERAGE = null;

    public CollectGridssMetrics() {
    	// By default, only run metrics required by the GRIDSS pre-processing step
    	PROGRAM = new LinkedHashSet<>(Arrays.asList(Program.CollectInsertSizeMetrics));
    }
    
    public static void main(final String[] args) {
        new CollectGridssMetrics().instanceMainWithExit(args);
    }
    @Override
    protected String[] customCommandLineValidation() {
        if (PROGRAM.isEmpty() && GRIDSS_PROGRAM.isEmpty()) {
            return new String[]{"No programs specified with PROGRAM or GRIDSS_PROGRAM"};
        }
        if (GRIDSS_PROGRAM.contains(GridssProgram.ReportThresholdCoverage) && THRESHOLD_COVERAGE == null) {
        	return new String[]{"THRESHOLD_COVERAGE is required when running ReportThresholdCoverage."};
        }
        return super.customCommandLineValidation();
    }
    public int doWork() {
    	List<ProgramInterface> toRun = Lists.newArrayList(Iterables.transform(GRIDSS_PROGRAM, p -> new GridssProgramProgramInterfaceFactory().create(p)));
		toRun.addAll(PROGRAM);
    	setProgramsToRun(toRun);
    	if (REFERENCE_SEQUENCE != null && REFERENCE_SEQUENCE.exists()) {
			try {
				ReferenceCommandLineProgram.ensureDictionaryMatches(INPUT, new FastaSequenceFile(REFERENCE_SEQUENCE, true).getSequenceDictionary(), REFERENCE_SEQUENCE);
			} catch (IOException e) {
				throw new RuntimeIOException(e);
			}
		}
    	return super.doWork();
    }
    protected class GridssProgramProgramInterfaceFactory {
    	public ProgramInterface create(GridssProgram program) {
    		switch (program) {
    			case CollectCigarMetrics:
    				return new ProgramInterface() {
	    				@Override
	    				public SinglePassSamProgram makeInstance(final String outbase,
																 final String outext,
																 final File input,
																 final File reference,
																 final Set<MetricAccumulationLevel> metricAccumulationLevel,
																 final File dbSnp,
																 final File intervals,
																 final File refflat,
																 final  Set<String> ignoreSequence) {
	    					final CollectCigarMetrics program = new CollectCigarMetrics();
							program.output = new RequiredOutputArgumentCollection(new File(outbase + gridss.analysis.CollectCigarMetrics.METRICS_SUFFIX + outext));
	
	    	                // Generally programs should not be accessing these directly but it might make things smoother
	    	                // to just set them anyway. These are set here to make sure that in case of a the derived class
	    	                // overrides
	    	                program.INPUT = input;
	    	                program.setReferenceSequence(reference);
	
	    	                return program;
	    				}
	    				@Override
	    				public boolean needsReferenceSequence() {
	    					return false;
	    				}
	    				@Override
	    				public boolean supportsMetricAccumulationLevel() {
	    					return false;
	    				}
	    	        };
    			case CollectIdsvMetrics:
    				return new ProgramInterface() {
    					@Override
    					public SinglePassSamProgram makeInstance(final String outbase,
																 final String outext,
																 final File input,
																 final File reference,
																 final Set<MetricAccumulationLevel> metricAccumulationLevel,
																 final File dbSnp,
																 final File intervals,
																 final File refflat,
																 final  Set<String> ignoreSequence) {
    						final CollectIdsvMetrics program = new CollectIdsvMetrics();
							program.output = new RequiredOutputArgumentCollection(new File(outbase + gridss.analysis.CollectIdsvMetrics.METRICS_SUFFIX + outext));

    		                // Generally programs should not be accessing these directly but it might make things smoother
    		                // to just set them anyway. These are set here to make sure that in case of a the derived class
    		                // overrides
    		                program.INPUT = input;
    		                program.setReferenceSequence(reference);
    		                return program;
    					}

    					@Override
    					public boolean needsReferenceSequence() {
    						return false;
    					}
    					@Override
    					public boolean supportsMetricAccumulationLevel() {
    						return false;
    					}
    		        };
    			case CollectMapqMetrics:
    				return new ProgramInterface() {
    					@Override
    					public SinglePassSamProgram makeInstance(final String outbase,
																 final String outext,
																 final File input,
																 final File reference,
																 final Set<MetricAccumulationLevel> metricAccumulationLevel,
																 final File dbSnp,
																 final File intervals,
																 final File refflat,
																 final  Set<String> ignoreSequence) {
    						final CollectMapqMetrics program = new CollectMapqMetrics();
							program.output = new RequiredOutputArgumentCollection(new File(outbase + gridss.analysis.CollectMapqMetrics.METRICS_SUFFIX + outext));
    		                // TODO: why does this only work on my local desktop+servers and not the CI or cluster systems?
    		                //program.Histogram_FILE = new File(outbase + gridss.analysis.CollectMapqMetrics.HISTOGRAM_SUFFIX);

    		                // Generally programs should not be accessing these directly but it might make things smoother
    		                // to just set them anyway. These are set here to make sure that in case of a the derived class
    		                // overrides
    		                program.INPUT = input;
    		                program.setReferenceSequence(reference);

    		                return program;
    					}

    					@Override
    					public boolean needsReferenceSequence() {
    						return false;
    					}

    					@Override
    					public boolean supportsMetricAccumulationLevel() {
    						return false;
    					}
    		        };
    			case CollectTagMetrics:
    				return new ProgramInterface() {
    					@Override
    					public SinglePassSamProgram makeInstance(final String outbase,
																 final String outext,
																 final File input,
																 final File reference,
																 final Set<MetricAccumulationLevel> metricAccumulationLevel,
																 final File dbSnp,
																 final File intervals,
																 final File refflat,
																 final  Set<String> ignoreSequence) {
    						final CollectTagMetrics program = new CollectTagMetrics();
							program.output = new RequiredOutputArgumentCollection(new File(outbase + gridss.analysis.CollectTagMetrics.METRICS_SUFFIX + outext));

    		                // Generally programs should not be accessing these directly but it might make things smoother
    		                // to just set them anyway. These are set here to make sure that in case of a the derived class
    		                // overrides
    		                program.INPUT = input;
    		                program.setReferenceSequence(reference);

    		                return program;
    					}

    					@Override
    					public boolean needsReferenceSequence() {
    						return false;
    					}

    					@Override
    					public boolean supportsMetricAccumulationLevel() {
    						return false;
    					}
    		        };
    			case ReportThresholdCoverage:
    				return new ProgramInterface() {
    					@Override
    					public SinglePassSamProgram makeInstance(final String outbase,
																 final String outext,
																 final File input,
																 final File reference,
																 final Set<MetricAccumulationLevel> metricAccumulationLevel,
																 final File dbSnp,
																 final File intervals,
																 final File refflat,
																 final  Set<String> ignoreSequence) {
    						final ReportThresholdCoverage program = new ReportThresholdCoverage();
							program.output = new RequiredOutputArgumentCollection(new File(outbase + gridss.analysis.ReportThresholdCoverage.SUFFIX + outext));
    		                
    		                // Generally programs should not be accessing these directly but it might make things smoother
    		                // to just set them anyway. These are set here to make sure that in case of a the derived class
    		                // overrides
    		                program.INPUT = input;
    		                program.setReferenceSequence(reference);
    		                program.THRESHOLD_COVERAGE = THRESHOLD_COVERAGE;
    		                return program;
    					}

    					@Override
    					public boolean needsReferenceSequence() {
    						return false;
    					}

    					@Override
    					public boolean supportsMetricAccumulationLevel() {
    						return false;
    					}
    		        };
    		}
    		throw new NotImplementedException("Unreachable code encountered");
    	}
    }
}
