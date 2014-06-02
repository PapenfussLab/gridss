package au.edu.wehi.idsv.sam;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.nio.charset.StandardCharsets;

import au.edu.wehi.idsv.BreakendDirection;

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
	public String getBreakpointBases() {
		return new String(getReadBases(), getDirection() == BreakendDirection.Backward ? 0 : getAnchorLength(), getBreakpointLength(), StandardCharsets.US_ASCII);
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
	public float getAssemblyBreakpointQuality() {
		byte[] qual = getBaseQualities();
		int length = getBreakpointLength();
		if (length == 0) return 0;
		int start = getDirection() == BreakendDirection.Backward ? 0 : getAnchorLength();  
		int qualSum = 0;
		for (int i = 0; i < length; i++) {
			qualSum += qual[start + i];
		}
		return qualSum/(float)length;
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
