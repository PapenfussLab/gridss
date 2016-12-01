package gridss;

import java.util.Iterator;
import java.util.stream.Collectors;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.GreedyVariantAllocationCache;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IdsvVariantContextBuilder;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;

@CommandLineProgramProperties(
        usage = "Evidence reallocation is required to ensure that any given read/read pair/assembly, "
        		+ " supports only a consistent set breakpoints. \n"
        		+ "For discordant read pairs, this means that only a single breakpoint can be supported. \n"
        		+ "For split reads, only a single set of partial read mappings can be supported, and only one breakpoint per split can be supported.\n"
        		+ "For indels, only a single read mapping can be supported, and only one breakpoint per indel can be supported.\n"
        		+ "Note: the EID attribute must be populated with the relevant GRIDSS EvidenceID's using AllocateEvidence",  
        usageShort = "Unique allocates evidence supporting multiple mutually exclusive breakpoints."
)
public class AllocateEvidence extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AllocateEvidence.class);
	public static enum EvidenceAllocationStrategy {
		GREEDY,
	}
	private GreedyVariantAllocationCache cache;
	@Option(doc="Evidence allocation strategy used to uniquely assign evidence.")
	public EvidenceAllocationStrategy ALLOCATION_STRATEGY = EvidenceAllocationStrategy.GREEDY;
	@Override
	public CloseableIterator<VariantContextDirectedBreakpoint> iterator(CloseableIterator<VariantContextDirectedBreakpoint> calls) {
		populateCache();
		Iterator<VariantContextDirectedBreakpoint> it = Iterators.transform(calls, bp -> annotate(bp));
		return new AutoClosingIterator<>(it, calls);
	}
	public void populateCache() {
		IOUtil.assertFileIsReadable(INPUT_VCF);
		try (CloseableIterator<VariantContextDirectedBreakpoint> bpit = getBreakpoints(INPUT_VCF)) {
			populateCache(bpit);
		}
	}
	public void populateCache(CloseableIterator<VariantContextDirectedBreakpoint> calls) {
		log.info("Caching breakpoint evidence allocation");
		cache = new GreedyVariantAllocationCache();
		while (calls.hasNext()) {
			VariantContextDirectedBreakpoint bp = calls.next();
			String variantId = bp.getID();
			for (Object evidenceID : bp.getAttributeAsList(VcfAttributes.EVIDENCE_ID.attribute())) {
				cache.addBreakpoint(variantId, bp.getBreakpointQual(), (String)evidenceID);
			}
		}
		log.info("Caching breakpoint evidence allocation complete");
	}
	private VariantContextDirectedBreakpoint annotate(VariantContextDirectedBreakpoint bp) {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext(), bp);
		String variantId = bp.getID();
		float variantQual = bp.getBreakpointQual();
		builder.attribute(VcfAttributes.EVIDENCE_ID,
				bp.getAttributeAsList(VcfAttributes.EVIDENCE_ID.attribute()).stream()
				.filter(evidenceID -> cache.isBestBreakpoint(variantId, variantQual, (String)evidenceID))
				.collect(Collectors.toList()));
		IdsvVariantContext variant = builder.make();
		return (VariantContextDirectedBreakpoint)variant;
	}
}
