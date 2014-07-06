package au.edu.wehi.idsv.sam;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.SAMRecordUtil;

/**
 * Breakend assembly
 * 
 * Assembly information is encoded in SAMRecord to enable easy visualisation and
 * debugging of the assembly. Note that the genomic position of the assembly is
 * not encoded in the SAMRecord.
 * @author Daniel Cameron
 *
 */
public class AnomolousReadAssembly extends SAMRecord {
	public AnomolousReadAssembly(String assemblerProgram, byte[] assemblyBases, int assemblyAnchorLength, BreakendDirection direction) {
		this(assemblerProgram, assemblyBases, null, assemblyAnchorLength, direction, null, null);
	}
	public AnomolousReadAssembly(String assemblerProgram, byte[] assemblyBases, byte[] assemblyBaseQuality, int assemblyAnchorLength, BreakendDirection direction, Integer assembledReadCount, Integer assembledReadBaseCount) {
		super(null);
		setAttribute(SamTags.ASSEMBLER_PROGRAM, assemblerProgram);
		if (assembledReadCount != null) setAttribute(SamTags.ASSEMBLER_READ_COUNT, assembledReadCount);
		if (assembledReadBaseCount != null) setAttribute(SamTags.ASSEMBLER_READ_BASE_COUNT, assembledReadBaseCount);
		setReadBases(assemblyBases);
		setBaseQualities(assemblyBaseQuality);
		int assemblyLength = assemblyBases.length;
		if (direction == BreakendDirection.Forward) {
			setCigarString(String.format("%dM%dS", assemblyAnchorLength, assemblyLength - assemblyAnchorLength));
		} else {
			setCigarString(String.format("%dS%dM", assemblyLength - assemblyAnchorLength, assemblyAnchorLength));
		}
	}
	public String getAssemblerProgram() {
		return getStringAttribute(SamTags.ASSEMBLER_PROGRAM);
	}
	public BreakendDirection getDirection() {
		return getCigar().getCigarElement(0).getOperator() == CigarOperator.S ? BreakendDirection.Backward : BreakendDirection.Forward; 
	}
	public byte[] getBreakpointBases() {
		return getDirection() == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBases(this) : SAMRecordUtil.getStartSoftClipBases(this);
	}
	public byte[] getBreakpointQualities() {
		return getDirection() == BreakendDirection.Forward ? SAMRecordUtil.getEndSoftClipBaseQualities(this) : SAMRecordUtil.getStartSoftClipBaseQualities(this);
	}
	public String getAnchorBases() {
		return new String(getReadBases(), getDirection() == BreakendDirection.Forward ? 0 : getBreakpointLength(), getAnchorLength(), StandardCharsets.US_ASCII);
	}
	public int getAnchorLength() {
		if (getDirection() == BreakendDirection.Forward) return getCigar().getCigarElement(0).getLength();
		else return getCigar().getCigarElement(1).getLength(); 
	}
	public int getBreakpointLength() {
		if (getDirection() == BreakendDirection.Forward) return getCigar().getCigarElement(1).getLength();
		else return getCigar().getCigarElement(0).getLength(); 
	}
	/**
	 * Get average base quality of assembled reads
	 * @return
	 */
	public double getAverageBreakpointQuality() {
		byte[] quals = getBreakpointQualities();
		if (quals == null || quals.length == 0) return 0;
		int qualSum = 0;
		for (byte qual : quals) {
			qualSum += qual;
		}
		return qualSum / (float)quals.length;
	}
	/**
	 * Number of reads in assembly
	 * @return
	 */
	public Integer getReadCount() {
		return (Integer)getAttribute(SamTags.ASSEMBLER_READ_COUNT);
	}
	public Integer getReadBaseCount() {
		return (Integer)getAttribute(SamTags.ASSEMBLER_READ_BASE_COUNT);
	}
}
