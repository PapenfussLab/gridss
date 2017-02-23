package au.edu.wehi.idsv;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Iterator;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * VCF Breakend record
 * see Section 5.4.9 of http://samtools.github.io/hts-specs/VCFv4.2.pdf for details of breakends
 * @author Daniel Cameron
 *
 */
public class VariantContextDirectedEvidence extends IdsvVariantContext implements DirectedEvidence {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final VcfBreakendSummary breakend;
	//private static Log LOG = Log.getInstance(VariantContextDirectedBreakpoint.class);
	public VariantContextDirectedEvidence(GenomicProcessingContext processContext, EvidenceSource source, VariantContext context) {
		super(processContext, source, context);
		this.breakend = new VcfBreakendSummary(processContext, context);
	}
	@Override
	public BreakendSummary getBreakendSummary() {
		if (breakend.location == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakend.location;
	}
	@Override
	public String getEvidenceID() {
		return getID();
	}
	@Override
	public byte[] getBreakendSequence() {
		return breakend.breakpointSequence.getBytes(StandardCharsets.US_ASCII);
	}
	public String getBreakpointSequenceString() {
		if (breakend.breakpointSequence == null) throw new IllegalStateException(String.format("%s not a valid breakend", getID()));
		return breakend.breakpointSequence;
	}
	@Override
	public byte[] getBreakendQuality() {
		return null;
	}
	@Override
	public byte[] getAnchorSequence() {
		return null;
	}
	@Override
	public byte[] getAnchorQuality() {
		return null;
	}
	@Override
	public boolean isValid() {
		return breakend.location != null;
	}
	/**
	 * Returns an iterator containing only the breakend variants from the given iterator
	 * @param context processing context
	 * @param variants input variants
	 * @return breakend variants
	 */
	public static Iterator<VariantContextDirectedEvidence> breakendIterator(final ProcessingContext context, final EvidenceSource source, final Iterator<VariantContext> variants) {
		 Iterator<IdsvVariantContext> itivc = Iterators.transform(variants, new Function<VariantContext, IdsvVariantContext>() {
			public IdsvVariantContext apply(VariantContext v) {
				return new IdsvVariantContextBuilder(context, v).source(source).make();
			}
		});
		 Iterator<VariantContextDirectedEvidence> itde = Iterators.filter(Iterators.filter(itivc, VariantContextDirectedEvidence.class), new Predicate<VariantContextDirectedEvidence>() {
			public boolean apply(VariantContextDirectedEvidence v) {
				return v.isValid();
			}
		});
		return itde;
	}
	public static Ordering<VariantContextDirectedEvidence> ByBreakendStartEnd = new Ordering<VariantContextDirectedEvidence>() {
		public int compare(VariantContextDirectedEvidence o1, VariantContextDirectedEvidence o2) {
			return BreakendSummary.ByStartEnd.compare(o1.getBreakendSummary(), o2.getBreakendSummary());
		  }
	};
	@Override
	public int getLocalMapq() {
		throw new IllegalArgumentException("NYI");
	}
	@Override
	public boolean isBreakendExact() {
		return !hasAttribute(VcfSvConstants.IMPRECISE_KEY);
	}
	@Override
	public float getBreakendQual() {
		return (float)getPhredScaledQual();
	}
	/**
	 * Converts to SAMRecord representation.
	 * @param header
	 * @return SAMRecord
	 */
	public SAMRecord asSamRecord(SAMFileHeader header) {
		BreakendSummary be = getBreakendSummary();
		SAMRecord r = new SAMRecord(header);
		Cigar cigar = new Cigar(be.getCigarRepresentation());
		if (hasAttribute(VcfInfoAttributes.SUPPORT_CIGAR.attribute())) {
			cigar = TextCigarCodec.decode(getAttributeAsString(VcfInfoAttributes.SUPPORT_CIGAR.attribute(), null));
		}
		r.setCigar(cigar);
		r.setReadName(getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, getID()));
		byte[] bases = new byte[cigar.getReadLength()];
		Arrays.fill(bases, (byte)'N');
		r.setReadBases(bases);
		r.setBaseQualities(SAMRecord.NULL_SEQUENCE);
		r.setReferenceIndex(be.referenceIndex);
		if (be.direction == BreakendDirection.Forward) {
			r.setReadNegativeStrandFlag(false);
			r.setAlignmentStart(be.end - cigar.getReadLength() + 1);
		} else {
			r.setReadNegativeStrandFlag(true);
			r.setAlignmentStart(getBreakendSummary().start);
		}
		if (be instanceof BreakpointSummary) {
			BreakendSummary rbe = ((BreakpointSummary)be).remoteBreakend();
			r.setReadPairedFlag(true);
			r.setMateNegativeStrandFlag(rbe.direction == BreakendDirection.Backward);
			r.setMateReferenceIndex(rbe.referenceIndex);
			// incorrect if mate is FWD and has an ANCHOR_CIGAR defined
			r.setMateAlignmentStart(rbe.start);
			r.setMateUnmappedFlag(false);
			r.setFirstOfPairFlag(getID().endsWith("o"));
			r.setSecondOfPairFlag(getID().endsWith("h"));
		}
		return r;
	}
	public int getBreakendEvidenceCountAssembly() { return getInt(VcfInfoAttributes.BREAKEND_ASSEMBLY_COUNT, 0); }
	public int getBreakendEvidenceCountReadPair() { return getInt(VcfInfoAttributes.BREAKEND_UNMAPPEDMATE_COUNT, 0); }
	public int getBreakendEvidenceCountReadPair(int category) { return getInt(category, VcfFormatAttributes.BREAKEND_UNMAPPEDMATE_COUNT, 0); } 
	public int getBreakendEvidenceCountSoftClip() { return getInt(VcfInfoAttributes.BREAKEND_SOFTCLIP_COUNT, 0); }
	public int getBreakendEvidenceCountSoftClip(int category) { return getInt(category, VcfFormatAttributes.BREAKEND_SOFTCLIP_COUNT, 0); }
	public int getReferenceReadCount() { return getInt(VcfInfoAttributes.REFERENCE_READ_COUNT, 0); }
	public int getReferenceReadCount(int category) { return getInt(category, VcfFormatAttributes.REFERENCE_READ_COUNT, 0); }
	public int getReferenceReadPairCount() { return getInt(VcfInfoAttributes.REFERENCE_READPAIR_COUNT, 0); }
	public int getReferenceReadPairCount(int category)  { return getInt(category, VcfFormatAttributes.REFERENCE_READPAIR_COUNT, 0); }
	/**
	 * Assembled breakend sequence. 
	 * @return Breakend assembly sequence, null if no assembly was found
	 */
	public FastqRecord getBreakendAssemblyFastq() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public boolean isFromMultimappingFragment() {
		throw new UnsupportedOperationException();
	}
}
