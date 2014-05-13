package au.edu.wehi.socrates;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Processing context for the given record
 * @author Daniel Cameron
 *
 */
public class ProcessingContext {
	private final ReferenceSequenceFile reference;
	private final SAMSequenceDictionary dictionary;
	private final LinearGenomicCoordinate linear;
	private final RelevantMetrics metrics;
	public ProcessingContext(
		ReferenceSequenceFile reference,
		SAMSequenceDictionary dictionary,
		RelevantMetrics metrics) {
		this.reference = reference;
		this.dictionary = dictionary;
		this.linear = new LinearGenomicCoordinate(dictionary);
		this.metrics = metrics;
	}
	public ProcessingContext(
			ReferenceSequenceFile reference,
			RelevantMetrics metrics) {
			this.reference = reference;
			this.dictionary = reference.getSequenceDictionary();
			if (this.dictionary == null) {
				throw new RuntimeException("Missing sequence dictionary for reference genome. Please create using picard CreateSequenceDictionary.");
			}
			this.linear = new LinearGenomicCoordinate(dictionary);
			this.metrics = metrics;
		}
	public ReferenceSequenceFile getReference() {
		return reference;
	}
	public SAMSequenceDictionary getDictionary() {
		return dictionary;
	}
	public LinearGenomicCoordinate getLinear() {
		return linear;
	}
	public RelevantMetrics getMetrics() {
		return metrics;
	}
}
