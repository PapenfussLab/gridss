package au.edu.wehi.socrates;

import java.util.Iterator;

import net.sf.samtools.SAMRecord;

public class SAMRecordEvidenceProcessorOrchestrator {
	public SAMRecordEvidenceProcessorOrchestrator(
		Iterator<SAMRecord> sc,
		Iterator<SAMRecord> oea,
		Iterator<SAMRecord> dp,
		Iterator<SAMRecord> mate) {
	}
	public void process() {
		// genomic position = first
		// while input reads
			// flush expired reads from working sets
			// add reads now in working set
			// process genomic position
			// advance genomic position to next callable position
		// process remaining callable positions
	}
	
	public interface MethodConsumer<T> {
		T process(SAMRecord r);
	}
	public interface IterableConsumer<T> extends Iterable<T> {
	}
	
	// Directed Breakpoint Outputs:
	public interface DirectedBreakpointCalculator extends Iterable<T> {
	}
	
	
	// Orchestrator
	// construct working set
	// ask calculators to process positions
	// optimisation: don't ask
	
	// calculators need to tell orchestrator:
	// what reads are included in their working set:
	// {(Type, min, max)}
	
	// each record
	// - chain through every dirbpcalc
	//
	
	// backward processing:
	// for each record sorted by end
}
