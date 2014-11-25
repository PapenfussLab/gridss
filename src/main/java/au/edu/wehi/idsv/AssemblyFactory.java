package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import au.edu.wehi.idsv.util.CollectionUtil;
import au.edu.wehi.idsv.vcf.VcfAttributes;

import com.google.common.base.Function;
import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public final class AssemblyFactory {
	private AssemblyFactory() { } 
	/**
	 * Creates an assembly 
	 * @param processContext context
	 * @param source assembly source
	 * @param direction direction of breakend
	 * @param evidence evidence supporting the assembly breakend
	 * @param anchorReferenceIndex contig of anchored bases 
	 * @param anchorBreakendPosition genomic position of anchored base closest breakend
	 * @param anchoredBaseCount number of anchored bases in assembly
	 * @param baseCalls assembly base sequence as per a positive strand read over the anchor
	 * @param baseQuals assembly base qualities
	 * @param normalBaseCount number of assembly bases contributed by normal evidence sources
	 * @param tumourBaseCount number of assembly bases contributed by tumour evidence sources
	 * @return assembly evidence for the given assembly
	 */
	public static SAMRecordAssemblyEvidence createAnchored(
			ProcessingContext processContext,
			AssemblyEvidenceSource source, BreakendDirection direction,
			Set<DirectedEvidence> evidence,
			int anchorReferenceIndex, int anchorBreakendPosition, int anchoredBaseCount,
			byte[] baseCalls, byte[] baseQuals,
			int normalBaseCount, int tumourBaseCount) {
		BreakendSummary breakend = new BreakendSummary(anchorReferenceIndex, direction, anchorBreakendPosition, anchorBreakendPosition);
		Map<VcfAttributes, int[]> attr = calculateAssemblyIntAttributes(evidence, normalBaseCount, tumourBaseCount);
		if (attr.get(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT)[0] + attr.get(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT)[1] == 0 &&
				attr.get(VcfAttributes.ASSEMBLY_READPAIR_COUNT)[0] + attr.get(VcfAttributes.ASSEMBLY_READPAIR_COUNT)[1] > 0) {
			throw new RuntimeException(String.format("Sanity check failure: anchored assembly at %s created from only unanchored evidence", new BreakendSummary(anchorReferenceIndex, direction, anchorBreakendPosition, anchorBreakendPosition)));
		}
		return new SAMRecordAssemblyEvidence(processContext.getBasicSamHeader(), breakend, source, anchoredBaseCount, baseCalls, baseQuals, attr);
	}
	/**
	 * Creates an assembly whose breakpoint cannot be exactly anchored to the reference  
	 * @param processContext context
	 * @param source assembly source
	 * @param direction direction of breakend
	 * @param evidence evidence supporting the assembly breakend
	 * @param baseCalls assembly base sequence as per a positive strand read into a putative anchor
	 * @param baseQuals assembly base qualities
	 * @param normalBaseCount number of assembly bases contributed by normal evidence sources
	 * @param tumourBaseCount number of assembly bases contributed by tumour evidence sources
	 * @return assembly evidence for the given assembly
	 */
	public static SAMRecordAssemblyEvidence createUnanchored(
			ProcessingContext processContext,
			AssemblyEvidenceSource source,
			Set<DirectedEvidence> evidence,
			byte[] baseCalls, byte[] baseQuals,
			int normalBaseCount, int tumourBaseCount) {
		BreakendSummary breakend = Models.calculateBreakend(Lists.newArrayList(evidence));
		Map<VcfAttributes, int[]> attr = calculateAssemblyIntAttributes(evidence, normalBaseCount, tumourBaseCount);
		return new SAMRecordAssemblyEvidence(processContext.getBasicSamHeader(), breakend, source, 0, baseCalls, baseQuals, attr);
	}
	protected static Map<VcfAttributes, int[]> calculateAssemblyIntAttributes(Set<DirectedEvidence> evidence, int normalBaseCount, int tumourBaseCount) {
		if (evidence == null) {
			evidence = Sets.newHashSet();
		}
		List<NonReferenceReadPair> rp = Lists.newArrayList(Iterables.filter(evidence, NonReferenceReadPair.class));
		List<SoftClipEvidence> sc = Lists.newArrayList(Iterables.filter(evidence, SoftClipEvidence.class));
		List<NonReferenceReadPair> rpNormal = Lists.newArrayList(Iterables.filter(rp, new Predicate<NonReferenceReadPair>() { public boolean apply(NonReferenceReadPair e) { return !((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		List<NonReferenceReadPair> rpTumour = Lists.newArrayList(Iterables.filter(rp, new Predicate<NonReferenceReadPair>() { public boolean apply(NonReferenceReadPair e) { return ((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		List<SoftClipEvidence> scNormal = Lists.newArrayList(Iterables.filter(sc, new Predicate<SoftClipEvidence>() { public boolean apply(SoftClipEvidence e) { return !((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		List<SoftClipEvidence> scTumour = Lists.newArrayList(Iterables.filter(sc, new Predicate<SoftClipEvidence>() { public boolean apply(SoftClipEvidence e) { return ((SAMEvidenceSource)e.getEvidenceSource()).isTumour(); } }) );
		
		HashMap<VcfAttributes, int[]> attributes = Maps.newHashMap();
		attributes.put(VcfAttributes.ASSEMBLY_BASE_COUNT, new int[] { normalBaseCount, tumourBaseCount });		
		attributes.put(VcfAttributes.ASSEMBLY_READPAIR_COUNT, new int[] { rpNormal.size(), rpTumour.size() } );
		attributes.put(VcfAttributes.ASSEMBLY_SOFTCLIP_COUNT, new int[] { scNormal.size(), scTumour.size()} );
		Function<SoftClipEvidence, Integer> fscLen = new Function<SoftClipEvidence, Integer>() { public Integer apply(SoftClipEvidence e) { return e.getSoftClipLength(); } };
		Function<NonReferenceReadPair, Integer> frpReadLength = new Function<NonReferenceReadPair, Integer>() { public Integer apply(NonReferenceReadPair e) { return e.getNonReferenceRead().getReadLength(); } };
		int scLenN = CollectionUtil.maxInt(scNormal, fscLen, 0);
		int scLenT = CollectionUtil.maxInt(scTumour, fscLen, 0);
		//int scLen = Math.max(scLenN, scLenT);
		attributes.put(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_TOTAL, new int[] { CollectionUtil.sumInt(scNormal, fscLen), CollectionUtil.sumInt(scTumour, fscLen) } );
		attributes.put(VcfAttributes.ASSEMBLY_SOFTCLIP_CLIPLENGTH_MAX, new int[] { scLenN, scLenT } );
		int rpReadLenN = CollectionUtil.maxInt(rpNormal, frpReadLength, 0);
		int rpReadLenT = CollectionUtil.maxInt(rpTumour, frpReadLength, 0);
		//int rpReadLen = Math.max(rpReadLenN, rpReadLenT);
		attributes.put(VcfAttributes.ASSEMBLY_READPAIR_LENGTH_MAX, new int[] { rpReadLenN, rpReadLenT} );
		
		int localMapq = 0;
		for (DirectedEvidence e : evidence) {
			localMapq = Math.max(localMapq, e.getLocalMapq());
		}
		attributes.put(VcfAttributes.ASSEMBLY_MAPQ_LOCAL_MAX, new int[] { localMapq });
		return attributes;
	}
	/**
	 * Updates the given assembly to incorporate the given realignment of the assembly breakend
	 * @param processContext
	 * @return
	 */
	public static SAMRecordAssemblyEvidence incorporateRealignment(ProcessingContext processContext, SAMRecordAssemblyEvidence assembly, SAMRecord realignment) {
		if (realignment == null) return assembly;
		SAMRecordAssemblyEvidence a = (SAMRecordAssemblyEvidence)assembly;
		if (realignment.getReadUnmappedFlag() || !processContext.getRealignmentParameters().realignmentPositionUnique(realignment)) {
			// Breakend did not align well enough for us to call a breakpoint
			return new SAMRecordAssemblyEvidence(a.getEvidenceSource(), a.getSAMRecord(), realignment);
		} else {
			return new RealignedSAMRecordAssemblyEvidence(processContext, a.getEvidenceSource(), a.getSAMRecord(), realignment);
		}
	}
}
