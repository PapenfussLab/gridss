package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.sam.ChimericAlignment;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;

public class GreedyVariantAllocationCache {
	/**
	 * read (pair) -> (event, score, read pair alignment)
	 * Only best placement of the read pairs should be allocated.
	 */
	private final HashMap<String, Node> bestReadPairAlignment = new HashMap<>();
	/**
	 * read -> (event, score, read alignment)
	 * 
	 * Only evidence from the best alignment for each read should be allocated
	 * since we don't want to allocate the mutually exclusive evidence that
	 * results from two separate read alignments for the same read
	 */
	private final HashMap<String, Node> bestReadAlignment = new HashMap<>();
	/**
	 * evidenceID -> best event lookup
	 * Each evidence can support only a single variant.
	 */
	private final HashMap<String, Node> bestEventForEvidence = new HashMap<>();
	public void addBreakpoint(String event, float score, DirectedEvidence evidence) {
		put(bestEventForEvidence, evidence.getEvidenceID(), null, event, score);
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			String readpairid = dp.getLocalledMappedRead().getReadName();
			String alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			put(bestReadPairAlignment, readpairid, alignment, event, score);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			String readid = r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r));
			String alignment = getReadAlignment(r);
			put(bestReadAlignment, readid, alignment, event, score);
		}
	}
	private String getReadAlignment(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return "";
		List<ChimericAlignment> ca = Lists.newArrayList(ChimericAlignment.getChimericAlignments(r));
		ca.add(new ChimericAlignment(r));
		return getAlignmentString(ca);
	}
	private String getReadPairAlignment(SAMRecord r) {
		List<ChimericAlignment> ca = new ArrayList<>(2);
		if (!r.getReadUnmappedFlag()) {
			ca.add(new ChimericAlignment(r));
		}
		if (!r.getMateUnmappedFlag()) {
			ca.add(new ChimericAlignment(
					r.getMateReferenceName(),
					r.getMateAlignmentStart(),
					r.getMateNegativeStrandFlag(),
					TextCigarCodec.decode(r.getStringAttribute(SAMTag.MC.name())),
					// not using these fields so zeros are fine
					0, 
					0));
		}
		return getAlignmentString(ca);
	}
	private static final Comparator<ChimericAlignment> ByPosition = Comparator.<ChimericAlignment, Integer>comparing(ca -> ca.pos)
			.thenComparing(ca -> ca.isNegativeStrand)
			.thenComparing(ca -> ca.rname)
			.thenComparing(ca -> ca.cigar.toString());
	private String getAlignmentString(List<ChimericAlignment> ca) {
		StringBuilder sb = new StringBuilder();
		// need to make sure the order is the same at every alignment
		Collections.sort(ca, ByPosition);
		for (ChimericAlignment a : ca) {
			sb.append(a.rname);
			sb.append('#');
			sb.append(a.isNegativeStrand ? '-' : '+');
			sb.append(a.pos);
			sb.append('#');
			sb.append(a.cigar.toString());
			sb.append('#');
		}
		return sb.toString();
	}
	private static String getEvent(VariantContextDirectedBreakpoint variant) {
		return variant.getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, null);
	}
	public void addBreakpoint(VariantContextDirectedBreakpoint variant, DirectedEvidence evidence) {
		addBreakpoint(getEvent(variant), variant.getBreakpointQual(), evidence);
	}
	private boolean put(HashMap<String, Node> bestAlignment, String readName, String alignment, String event, float score) {
		Node existing = bestAlignment.get(readName);
		if (existing == null || existing.score < score) {
			bestAlignment.put(readName, new Node(event, score, alignment));
			return true;
		}
		return false;
	}
	private boolean isBestAlignment(HashMap<String, Node> bestAlignment, String readName, String alignment) {
		Node node = bestAlignment.get(readName);
		return node != null && Objects.equals(alignment, node.alignment);
	}
	public boolean isBestBreakpoint(String event, DirectedEvidence evidence) {
		if (!event.equals(bestEventForEvidence.get(evidence.getEvidenceID()).event)) {
			// This is not the best breakpoint supported by this evidence
			return false;
		}
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			String readpairid = dp.getLocalledMappedRead().getReadName();
			String alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			return isBestAlignment(bestReadPairAlignment, readpairid, alignment);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			String readid = r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r));
			String alignment = getReadAlignment(r);
			return isBestAlignment(bestReadAlignment, readid, alignment);
		}
	}
	public boolean isBestBreakpoint(VariantContextDirectedBreakpoint variant, DirectedEvidence evidence) {
		return isBestBreakpoint(getEvent(variant), evidence);
	}
	private static class Node {
		public Node(String event, float score, String alignment) {
			this.event = event;
			this.score = score;
			this.alignment = alignment;
		}
		public final String event;
		public final float score;
		public final String alignment;
		public String toString() {
			return String.format("(%s,%f,%s)", event, score, alignment);
		}
	}
}
