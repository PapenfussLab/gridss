package au.edu.wehi.idsv.sim;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;

import au.edu.wehi.idsv.ProcessingContext;

public class SimulatedChromosome {
	protected final byte[] seq;
	protected final int referenceIndex;
	protected final int margin;
	protected final ProcessingContext context;
	/**
	 * @param reference reference genome
	 * @param breakCleanMargin number of unambiguous bases around the breakpoint
	 */
	public SimulatedChromosome(ProcessingContext context, String chr, int margin) {
		this.context = context;
		this.referenceIndex = context.getReference().getSequenceDictionary().getSequence(chr).getSequenceIndex();
		this.seq = context.getReference().getSequence(chr).getBases();
		this.margin = margin;
	}
	protected void writeVcf(File vcf, Iterable<VariantContext> calls) {
		VariantContextWriter writer = context.getVariantContextWriter(vcf, true);
		for (VariantContext vc : calls) {
			writer.add(vc);
		}
		writer.close();
	}
	protected String getChr() { return context.getReference().getSequenceDictionary().getSequence(referenceIndex).getSequenceName(); }
}