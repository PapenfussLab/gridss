package gridss.analysis;

import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;

import picard.analysis.CollectMultipleMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

/**
 * Class that is designed to instantiate and execute multiple metrics programs that extend
 * SinglePassSamProgram while making only a single pass through the SAM file and supplying
 * each program with the records as it goes.
 *
 */
@CommandLineProgramProperties(
        usage = "Extension of picard.CollectMultipleMetrics to include gridss metrics. "
        		+ "Additional metrics are CollectCigarMetrics, CollectIdsvMetrics, and CollectMapqMetrics.",
        usageShort = "A \"meta-metrics\" calculating program that produces multiple metrics for the provided SAM/BAM",
        programGroup = Metrics.class
)
public class CollectGridssMetrics extends CollectMultipleMetrics {
    public static enum GridssProgram implements ProgramInterface {
    	CollectCigarMetrics {
			@Override
			public SinglePassSamProgram makeInstance(String outbase, String outext, File input, File reference,
					Set<MetricAccumulationLevel> metricAccumulationLevel, File dbSnp, File intervals) {
				final CollectCigarMetrics program = new CollectCigarMetrics();
                program.OUTPUT = new File(outbase + gridss.analysis.CollectCigarMetrics.METRICS_SUFFIX + outext);

                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

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
        },
    	CollectMapqMetrics {
			@Override
			public SinglePassSamProgram makeInstance(String outbase, String outext, File input, File reference,
					Set<MetricAccumulationLevel> metricAccumulationLevel, File dbSnp, File intervals) {
				final CollectMapqMetrics program = new CollectMapqMetrics();
                program.OUTPUT = new File(outbase + gridss.analysis.CollectMapqMetrics.METRICS_SUFFIX + outext);
                program.Histogram_FILE = new File(outbase + gridss.analysis.CollectMapqMetrics.HISTOGRAM_SUFFIX);

                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

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
        },
    	CollectTagMetrics {
			@Override
			public SinglePassSamProgram makeInstance(String outbase, String outext, File input, File reference,
					Set<MetricAccumulationLevel> metricAccumulationLevel, File dbSnp, File intervals) {
				final CollectTagMetrics program = new CollectTagMetrics();
                program.OUTPUT = new File(outbase + gridss.analysis.CollectTagMetrics.METRICS_SUFFIX + outext);

                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

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
        },
        CollectIdsvMetrics {
			@Override
			public SinglePassSamProgram makeInstance(String outbase, String outext, File input, File reference,
					Set<MetricAccumulationLevel> metricAccumulationLevel, File dbSnp, File intervals) {
				final CollectIdsvMetrics program = new CollectIdsvMetrics();
                program.OUTPUT = new File(outbase + gridss.analysis.CollectIdsvMetrics.METRICS_SUFFIX + outext);

                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.REFERENCE_SEQUENCE = reference;

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
        },
    }
    @Option(doc = "Set of gridss metrics programs to apply during the pass through the SAM file.")
    public Set<GridssProgram> GRIDSS_PROGRAM = new LinkedHashSet<>(Arrays.asList(GridssProgram.values()));

    public CollectGridssMetrics() {
    	// By default, only run those required by GRIDSS
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
        return super.customCommandLineValidation();
    }
    @Override
    public int doWork() {
    	List<ProgramInterface> toRun = Lists.newArrayList(GRIDSS_PROGRAM);
		toRun.addAll(PROGRAM);
    	setProgramsToRun(toRun);
    	return super.doWork();
    }
}
