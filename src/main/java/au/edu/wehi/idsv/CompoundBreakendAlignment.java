package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;

import com.google.common.collect.ImmutableList;

/**
 * Alignment of a sequence contig that spans multiple segments
 * @author Daniel Cameron
 *
 */
public class CompoundBreakendAlignment {
	private final ProcessingContext processContext;
	private final NavigableMap<Integer, SAMRecord> mappedStartOffset = new TreeMap<Integer, SAMRecord>();
	/**
	 * Placeholder anchor SAMRecord containing the minimal information required about the breakpoint
	 */
	private final SAMRecord firstanchor;
	public static int getAnchorLength(SAMRecord anchor) {
		if (getBreakendDirection(anchor) == BreakendDirection.Forward) {
			return anchor.getReadLength() - SAMRecordUtil.getEndSoftClipLength(anchor);
		} else {
			return anchor.getReadLength() - SAMRecordUtil.getStartSoftClipLength(anchor);
		}
	}
	public int getBreakpointCount() {
		return mappedStartOffset.size();
	}
	public SAMRecord getSimpleBreakendRealignment() {
		SAMRecord r;
		if (mappedStartOffset.size() == 0) {
			// no realignments
			r = getPlaceholderRealignment();
		} else {
			if (getBreakendDirection(firstanchor) == BreakendDirection.Forward) {
				r = mappedStartOffset.firstEntry().getValue();
			} else {
				r = mappedStartOffset.lastEntry().getValue();
			}
			r = asFullSequence(r);
		}
		trimRealignment(firstanchor, r);
		return r;
	}
	/**
	 * Determines where, relative to the realign sequence the anchored bases are 
	 * @param anchor anchoring
	 * @param realign realignment
	 * @return true if the anchored bases corresponding to the starting realign bases, false if they are at the end 
	 */
	private static boolean anchoredBasesAtStart(SAMRecord anchor, SAMRecord realign) {
		BreakendDirection anchorDirection = getBreakendDirection(anchor);
		boolean anchorForward = anchorDirection == BreakendDirection.Forward;
		boolean realignNegative = realign.getReadNegativeStrandFlag();
		if ((!anchorForward && !realignNegative) || (anchorForward && realignNegative)) {
			// BWD -> anchor comes from end
			// FWD + -ve
			return false;
		} else {
			// FWD -> anchor comes from start
			// BWD + -ve
			return true;
		}
	}
	private static BreakendDirection getBreakendDirection(SAMRecord anchor) {
		assert(!anchor.getReadNegativeStrandFlag()); // anchoring assemblies always considered to be on the +ve strand 
		Object anchorDirectionObject = anchor.getAttribute(SamTags.ASSEMBLY_DIRECTION);
		if (anchorDirectionObject instanceof Character) {
			return BreakendDirection.fromChar((char)(Character)anchorDirectionObject);
		} else {
			return null;
		}
	}
	/**
	 * Trims the given realignment based on the given anchor 
	 * @return
	 */
	private static void trimRealignment(SAMRecord anchor, SAMRecord realign) {
		if (anchoredBasesAtStart(anchor, realign)) {
			assert(realign.getReadUnmappedFlag() || SAMRecordUtil.getStartSoftClipLength(realign) >= getAnchorLength(anchor));
			SAMRecordUtil.trim(realign, getAnchorLength(anchor), 0);
		} else {
			assert(realign.getReadUnmappedFlag() || SAMRecordUtil.getEndSoftClipLength(realign) >= getAnchorLength(anchor));
			SAMRecordUtil.trim(realign, 0, getAnchorLength(anchor));
		}
	}
	/**
	 * Gets the set of subsequent breakpoints that the breakend sequence spans
	 * @return pair of SAMRecords with the left being the anchored alignment and the right the realigned 
	 */
	public List<Pair<SAMRecord, SAMRecord>> getSubsequentBreakpointAlignmentPairs() {
		if (getBreakpointCount() <= 1) return ImmutableList.of();
		List<Pair<SAMRecord, SAMRecord>> output = new ArrayList<Pair<SAMRecord,SAMRecord>>(getBreakpointCount() - 1);
		SAMRecord lastEntry = null;
		for (SAMRecord entry : getBreakendDirection(firstanchor) == BreakendDirection.Forward ? mappedStartOffset.values() : mappedStartOffset.descendingMap().values()) {
			if (lastEntry != null) {
				SAMRecord anchor = asFullSequence(lastEntry);
				SAMRecord realign = asFullSequence(entry);
				if (anchor.getReadNegativeStrandFlag()) {
					anchor.setReadNegativeStrandFlag(false); // assemblies are always considered to be +ve strand
					anchor.setAttribute(SamTags.ASSEMBLY_DIRECTION, getBreakendDirection(firstanchor).reverse().toChar());
					// need to flip realignment as well since we're aligning to the reverse comp
					realign.setReadNegativeStrandFlag(!realign.getReadNegativeStrandFlag());
				} else {
					anchor.setAttribute(SamTags.ASSEMBLY_DIRECTION, getBreakendDirection(firstanchor).toChar());
				}
				trimRealignment(anchor, realign);
				copyAssemblyAnnotations(anchor);
				anchor.setReadName(firstanchor.getReadName() + "_" + Integer.toString(output.size()));
				realign.setReadName(firstanchor.getReadName() + "_" + Integer.toString(output.size()));
				anchor.setMappingQuality(Math.min(anchor.getMappingQuality(), firstanchor.getMappingQuality()));
				output.add(Pair.of(anchor, realign));
			}
			lastEntry = entry;
		}
		return output;
	}
	private void copyAssemblyAnnotations(SAMRecord anchor) {
		for (String field: SamTags.ASSEMBLY_ANNOTATIONS) {
			if (field.equals(SamTags.ASSEMBLY_DIRECTION)) continue;
			anchor.setAttribute(field, firstanchor.getAttribute(field));
		}
	}
	private SAMRecord asFullSequence(SAMRecord realign) {
		SAMRecord r = SAMRecordUtil.clone(realign);
		assert(r.getReadLength() == r.getCigar().getReadLength()); // sanity check cigar matches read base count
		int readBreakendOffset = getReadOffset(realign);
		int startClipLength = SAMRecordUtil.getStartSoftClipLength(r);
		int endClipLength = SAMRecordUtil.getEndSoftClipLength(r);
		int length = r.getReadLength() - startClipLength - endClipLength;
		int startPad;
		int endPad;
		byte[] seq = firstanchor.getReadBases();
		byte[] qual = firstanchor.getBaseQualities();
		if (r.getReadNegativeStrandFlag()) {
			seq = Arrays.copyOf(seq, seq.length);
			qual = Arrays.copyOf(qual, qual.length);
			SequenceUtil.reverseComplement(seq);
			ArrayUtils.reverse(qual);
		}
		// A = anchor bases
		// # = start offset bases
		// S = start soft clip
		// E = end soft clip
		// ? = base counts that must be inferred
		// * = bases being aligned
		if (getBreakendDirection(firstanchor) == BreakendDirection.Forward) {
			// AAA###*********???
			if (r.getReadNegativeStrandFlag()) {
				// ???SSSMMMEEE###AAA
				endPad = endClipLength + readBreakendOffset + getAnchorLength(firstanchor);
				startPad = seq.length - endPad - length;
			} else {
				// AAAA###SSSMMMEEE???
				startPad = getAnchorLength(firstanchor) + readBreakendOffset + startClipLength;
				endPad = seq.length - startPad - length;
			}
		} else {
			// ###*********???AAA
			if (r.getReadNegativeStrandFlag()) {
				// AAA???SSSMMMEEE###
				endPad = endClipLength + readBreakendOffset;
				startPad = seq.length - endPad - length;
			} else {
				// ###SSSMMMEEE???AAA
				startPad = readBreakendOffset + startClipLength;
				endPad = seq.length - startPad - length;
			}
		}
		assert(startPad >= 0);
		assert(endPad >= 0);
		List<CigarElement> anchorCigarList = new ArrayList<CigarElement>(r.getCigar().getCigarElements());
		CigarUtil.trimClipping(anchorCigarList);
		CigarUtil.addStartSoftClip(anchorCigarList, startPad);
		CigarUtil.addEndSoftClip(anchorCigarList, endPad);
		r.setReadBases(seq);
		r.setBaseQualities(qual);
		r.setCigar(new Cigar(anchorCigarList));
		assert(r.getCigar().getReadLength() == seq.length);
		return r;
	}
	private SAMRecord getPlaceholderRealignment() {
		SAMRecord placeholder = new SAMRecord(firstanchor.getHeader());
		placeholder.setReadUnmappedFlag(true);
		placeholder.setReadBases(firstanchor.getReadBases());
		placeholder.setBaseQualities(firstanchor.getBaseQualities());
		placeholder.setReadNegativeStrandFlag(false);
		placeholder.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
		return placeholder;
	}
	private boolean isFiltered(SAMRecord realignment) {
		return realignment == null ||
				realignment.getReadUnmappedFlag() ||
				realignment.getMappingQuality() < processContext.getRealignmentParameters().mapqUniqueThreshold;
	}
	public CompoundBreakendAlignment(
			ProcessingContext processContext,
			SAMFileHeader header,
			BreakendSummary local,
			byte[] anchorSequence,
			byte[] anchorQual,
			byte[] breakendSequence,
			byte[] breakendQual,
			Collection<SAMRecord> realignments) {
		this(processContext, createPlaceholderAnchor(header, local, anchorSequence, anchorQual, breakendSequence, breakendQual), realignments);
	}
	public CompoundBreakendAlignment(
			ProcessingContext processContext,
			SAMRecord anchor,
			Collection<SAMRecord> realignments) {
		this.processContext = processContext;
		this.firstanchor = anchor;
		if (realignments != null) {
			for (SAMRecord r : realignments) {
				if (r != null && !isFiltered(r)) {
					mappedStartOffset.put(getMappedOffset(r), r);
				}
			}
		}
	}
	private static SAMRecord createPlaceholderAnchor(
			SAMFileHeader header,
			BreakendSummary local,
			byte[] anchorSequence,
			byte[] anchorQual,
			byte[] breakendSequence,
			byte[] breakendQual) {
		byte[] fullSeq;
		byte[] fullQual;
		if (local.direction == BreakendDirection.Forward) {
			fullSeq = ArrayUtils.addAll(anchorSequence, breakendSequence);
			fullQual = ArrayUtils.addAll(anchorQual, breakendQual);
		} else {
			fullSeq = ArrayUtils.addAll(breakendSequence, anchorSequence);
			fullQual = ArrayUtils.addAll(breakendQual, anchorQual);
		}
		SAMRecord r = new SAMRecord(header); 
		r.setReadBases(fullSeq);
		r.setBaseQualities(fullQual);
		//r.setReferenceIndex(local.referenceIndex);
		if (local.direction == BreakendDirection.Forward) {
			//r.setAlignmentStart(local.start);
			r.setCigar(new Cigar(ImmutableList.of(new CigarElement(anchorSequence.length, CigarOperator.M), new CigarElement(breakendSequence.length, CigarOperator.S))));
		} else {
			//r.setAlignmentStart(local.start - anchorSequence.length + 1);
			r.setCigar(new Cigar(ImmutableList.of(new CigarElement(breakendSequence.length, CigarOperator.S), new CigarElement(anchorSequence.length, CigarOperator.M))));
		}
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, local.direction.toChar());
		return r;
	}
	/**
	 * Gets the breakend offset of the first base
	 * @param r realignment
	 * @return breakend offset of the first read base
	 */
	private int getReadOffset(SAMRecord r) {
		return BreakpointFastqEncoding.getEncodedBreakendOffset(r.getReadName());
	}
	/**
	 * Gets the breakend offset of the first mapped base
	 * @param r realignment
	 * @return breakend offset of the first mapped read base
	 */
	private int getMappedOffset(SAMRecord r) {
		int offset = getReadOffset(r);
		if (r.getReadNegativeStrandFlag()) {
			offset += SAMRecordUtil.getEndSoftClipLength(r);
		} else {
			offset += SAMRecordUtil.getStartSoftClipLength(r);
		}
		return offset;
	}
}
