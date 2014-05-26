package au.edu.wehi.socrates;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.reference.ReferenceSequenceFile;

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
	private boolean vcf41mode = false;
	public ProcessingContext(
		ReferenceSequenceFile reference,
		SAMSequenceDictionary dictionary,
		RelevantMetrics metrics) {
		this.reference = reference;
		this.dictionary = new DynamicSAMSequenceDictionary(dictionary);
		this.linear = new LinearGenomicCoordinate(dictionary);
		this.metrics = metrics;
	}
	public ProcessingContext(
			ReferenceSequenceFile reference,
			RelevantMetrics metrics) {
			this.reference = reference;
			if (reference.getSequenceDictionary() == null) {
				throw new RuntimeException("Missing sequence dictionary for reference genome. Please create using picard CreateSequenceDictionary.");
			}
			this.dictionary = new DynamicSAMSequenceDictionary(reference.getSequenceDictionary());
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
	/**
	 * Determines whether VCF records should be compatible with VCF v4.1
	 */
	public boolean getVcf41Mode() {
		return vcf41mode;
	}
	public void setVcf41Mode(boolean compatable) {
		vcf41mode = compatable;
	}
}
