package au.edu.wehi.idsv;

import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.base.Function;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Calls simple variants from breakends
 * @author Daniel Cameron
 *
 */
public class BreakendToSimpleCallImpl {
	private static final Log log = Log.getInstance(BreakendToSimpleCallImpl.class);
	private final ProcessingContext processContext;
	private final int margin;
	private List<IdsvVariantContext> outputBuffer;
	/**
	 * Dirty hack so we can store VariantContextDirectedBreakpoint and use BreakendSummary objects to perform lookup 
	 */
	private static Ordering<Object> ByStart = new Ordering<Object>() {
		@Override
		public int compare(Object left, Object right) {
			BreakendSummary bsl;
			BreakendSummary bsr;
			String idr;
			String idl;
			if (left instanceof VariantContextDirectedBreakpoint) {
				bsl = ((VariantContextDirectedBreakpoint)left).getBreakendSummary();
			} else {
				bsl = (BreakendSummary)left;
			}
			if (right instanceof VariantContextDirectedBreakpoint) {
				bsr = ((VariantContextDirectedBreakpoint)right).getBreakendSummary();
			} else {
				bsr = (BreakendSummary)right;
			}
			if (left instanceof VariantContext) {
				idl = ((VariantContext)left).getID();
			} else {
				idl = left.toString();
			}
			if (right instanceof VariantContext) {
				idr = ((VariantContext)right).getID();
			} else {
				idr = right.toString();
			}
			return ComparisonChain.start()
			        .compare(bsl.referenceIndex, bsr.referenceIndex)
			        .compare(bsl.start, bsr.start)
			        .compare(idl, idr)
			        .result();
		}
	};
	private TreeSet<Object/*VariantContextDirectedBreakpoint*/> lookup;
	private HashMap<String, VariantContextDirectedBreakpoint> id;
	private TreeSet<VariantContextDirectedBreakpoint> byQual;
	private int maxWidth;
	public BreakendToSimpleCallImpl(ProcessingContext processContext) {
		this.processContext = processContext;
		this.margin = processContext.getVariantCallingParameters().breakendMargin;
	}
	public void convert(File breakendCalls, File simpleOutput) {
		try {
			maxWidth = 0;
			id = new HashMap<String, VariantContextDirectedBreakpoint>();
			outputBuffer = new ArrayList<IdsvVariantContext>();
			load(breakendCalls);
			process();
			write(simpleOutput);
		} catch (Exception e) {
			log.error(e);
		}
	}
	private void process() {
		byQual = new TreeSet<VariantContextDirectedBreakpoint>(IdsvVariantContext.ByQual.reverse().thenComparing(IdsvVariantContext.ByID));
		lookup = new TreeSet<Object>(ByStart);
		for (VariantContextDirectedBreakpoint bp : id.values()) {
			lookup.add(bp);
			if (bp.getID().endsWith("o")) {
				byQual.add(bp);
			} else {
				assert(bp.getID().endsWith("h"));
			}
		}
		while (!byQual.isEmpty()) {
			VariantContextDirectedBreakpoint bp = byQual.iterator().next(); 
			process(bp);
		}
	}
	private Collection<VariantContextDirectedBreakpoint> findOverlaps(BreakpointSummary be) {
		ArrayList<VariantContextDirectedBreakpoint> overlapping = new ArrayList<VariantContextDirectedBreakpoint>();
		for (Object obj : lookup.subSet(new BreakendSummary(be.referenceIndex, be.direction, be.start - maxWidth - margin - 1, be.start - maxWidth - margin - 1, Integer.MAX_VALUE),
				new BreakendSummary(be.referenceIndex, be.direction, Integer.MAX_VALUE, be.start + 2 * maxWidth + margin + 1, Integer.MAX_VALUE))) {
			VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)obj;
			if (processContext.getVariantCallingParameters().withMargin(bp.getBreakendSummary())
					.overlaps(processContext.getVariantCallingParameters().withMargin(be))) {
				overlapping.add(bp);
			}
		}
		Collections.sort(overlapping, IdsvVariantContext.ByQual.reverse());
		return overlapping;
	}
	/**
	 * Process the lower of a breakpoint pair
	 * @param bp
	 */
	private void process(VariantContextDirectedBreakpoint bp) {
		assert(bp.getID().endsWith("o"));
		remove(bp);
		VariantContextDirectedBreakpoint mate = id.get(bp.getAttribute(VcfSvConstants.MATE_BREAKEND_ID_KEY));
		if (mate == null) {
			log.warn(String.format("Breakpoint %s is missing mate", bp.getID()));
		}
		if (bp.getBreakendSummary().referenceIndex != bp.getBreakendSummary().referenceIndex2) {
			outputBuffer.add(bp);
			if (mate != null) {
				outputBuffer.add(mate);
			}
			return;
		}
		BreakpointSummary bs = bp.getBreakendSummary();
		// assert(bs.start <= bs.start2);
		// inversion
		if (bs.start < bs.start2 && bs.direction == bs.direction2) {
			VariantContextDirectedBreakpoint partner = null;
			for (VariantContextDirectedBreakpoint pairing : findOverlaps(new BreakpointSummary(
					bs.referenceIndex, bs.direction.reverse(), bs.start + 1, bs.start + 1, bs.end + 1, 
					bs.referenceIndex2, bs.direction2.reverse(), bs.start2 + 1, bs.start2 + 1, bs.end2 + 1))) {
				if (pairing != bp && pairing != mate) {
					partner = pairing;
					break;
				}
			}
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, bp);
			builder
				.phredScore(bp.getPhredScaledQual() + (partner == null ? 0 : partner.getPhredScaledQual()))
				.alleles("N", "<INV>")
				.start(bp.getBreakendSummary().start + 1)
				.stop(bs.start2)
				.attribute(VCFConstants.END_KEY, bs.start2)
				.attribute(VcfSvConstants.SV_TYPE_KEY, "INV");
				//.attribute(VcfSvConstants.SV_LENGTH_KEY, bs.start - bs.start2);
			setAttributes(bp, builder);
			IdsvVariantContext v = builder.make();
			assert(v != null);
			outputBuffer.add(v);
			remove(bp);
			if (partner != null) remove(partner);
			return;
		}
		if (bs.start < bs.start2 && bs.direction == BreakendDirection.Backward && bs.direction2 == BreakendDirection.Forward) {
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, bp);
			builder.alleles("N", "<DUP>")
				.start(bp.getBreakendSummary().start)
				.stop(bs.start2)
				.attribute(VCFConstants.END_KEY, bs.start2)
				.attribute(VcfSvConstants.SV_TYPE_KEY, "DUP");
				//.attribute(VcfSvConstants.SV_LENGTH_KEY, bs.start2 - bs.start);
			setAttributes(bp, builder);
			IdsvVariantContext v = builder.make();
			assert(v != null);
			outputBuffer.add(v);
			remove(bp);
			return;
		}
		if (bs.start < bs.start2 && bs.direction == BreakendDirection.Forward && bs.direction2 == BreakendDirection.Backward) {
			IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(processContext, bp);
			builder.alleles("N", "<DEL>")
				.start(bp.getBreakendSummary().start)
				.stop(bs.start2)
				.attribute(VCFConstants.END_KEY, bs.start2)
				.attribute(VcfSvConstants.SV_TYPE_KEY, "DEL");
				//.attribute(VcfSvConstants.SV_LENGTH_KEY, bs.start - bs.start2);
			setAttributes(bp, builder);
			IdsvVariantContext v = builder.make();
			assert(v != null);
			outputBuffer.add(v);
			remove(bp);
			return;
		}
		// just write out breakend
		outputBuffer.add(bp);
		if (mate != null) {
			outputBuffer.add(mate);
		}
		remove(bp);
	}
	private IdsvVariantContextBuilder setAttributes(VariantContextDirectedBreakpoint bp, IdsvVariantContextBuilder builder) {
		Object cirpos = bp.getAttribute(VcfInfoAttributes.CONFIDENCE_INTERVAL_REMOTE_BREAKEND_START_POSITION_KEY.attribute());
		if (cirpos != null) {
			builder.attribute(VcfSvConstants.CONFIDENCE_INTERVAL_END_POSITION_KEY, cirpos);
		}
		return builder;
	}
	private void remove(VariantContextDirectedBreakpoint bp) {
		lookup.remove(bp);
		byQual.remove(bp);
		VariantContextDirectedBreakpoint mate = id.get(bp.getAttribute(VcfSvConstants.MATE_BREAKEND_ID_KEY));
		if (mate != null) {
			lookup.remove(mate);
			byQual.remove(mate);
		} else {
			log.debug(String.format("%s missing mate", bp.getID()));
		}
	}
	private void write(File output) throws IOException {
		log.info("Writing simplified variant calls to " + output.getName() + ". Only use this output if your pipeline is unable to process variants in VCF breakend notation.");
		File working = FileSystemContext.getWorkingFileFor(output);
		VariantContextWriter vcfWriter = processContext.getVariantContextWriter(working, true);
		Collections.sort(outputBuffer, IdsvVariantContext.ByLocationStart);
		for (IdsvVariantContext v : outputBuffer) {
			vcfWriter.add(v);
		}
		CloserUtil.close(vcfWriter);
		FileHelper.move(working, output, true);
	}
	private void load(File breakendCalls) {
		log.info("Loading variants from " + breakendCalls.getName());
		int ignored = 0;
		VCFFileReader vcfReader = new VCFFileReader(breakendCalls, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		Iterator<IdsvVariantContext> idsvIt = Iterators.transform(it, new Function<VariantContext, IdsvVariantContext>() {
			@Override
			public IdsvVariantContext apply(VariantContext arg) {
				return IdsvVariantContext.create(processContext.getDictionary(), null, arg);
			}
		});
		while (idsvIt.hasNext()) {
			IdsvVariantContext v = idsvIt.next();
			if (v instanceof VariantContextDirectedBreakpoint) {
				VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint) v;
				id.put(v.getID(), bp);
				maxWidth = Math.max(maxWidth, bp.getBreakendSummary().end - bp.getBreakendSummary().start);
			} else {
				ignored++;
				assert(v != null);
				outputBuffer.add(v);
			}
		}
		it.close();
		vcfReader.close();
		if (ignored > 0) {
			log.warn(String.format("Ignored %d variants", ignored));
		}
	}
}
