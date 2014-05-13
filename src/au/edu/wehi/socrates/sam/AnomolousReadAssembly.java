package au.edu.wehi.socrates.sam;

import java.nio.charset.StandardCharsets;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import au.edu.wehi.socrates.BreakpointDirection;

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
	public AnomolousReadAssembly(String assemblerProgram, byte[] assemblyBases, int assemblyAnchorLength, BreakpointDirection direction) {
		this(assemblerProgram, assemblyBases, null, assemblyAnchorLength, direction, null);
	}
	public AnomolousReadAssembly(String assemblerProgram, byte[] assemblyBases, byte[] assemblyBaseQuality, int assemblyAnchorLength, BreakpointDirection direction, Integer assembledReadCount) {
		super(null);
		setAttribute(SocratesSamTags.ASSEMBLER_PROGRAM, assemblerProgram);
		if (assembledReadCount != null) setAttribute(SocratesSamTags.ASSEMBLER_READ_COUNT, assembledReadCount);
		setReadBases(assemblyBases);
		setBaseQualities(assemblyBaseQuality);
		int assemblyLength = assemblyBases.length;
		if (direction == BreakpointDirection.Forward) {
			setCigarString(String.format("%dM%dS", assemblyAnchorLength, assemblyLength - assemblyAnchorLength));
		} else {
			setCigarString(String.format("%dS%dM", assemblyLength - assemblyAnchorLength, assemblyAnchorLength));
		}
	}
	public String getAssemblerProgram() {
		return getStringAttribute(SocratesSamTags.ASSEMBLER_PROGRAM);
	}
	public BreakpointDirection getDirection() {
		return getCigar().getCigarElement(0).getOperator() == CigarOperator.S ? BreakpointDirection.Backward : BreakpointDirection.Forward; 
	}
	public String getBreakpointBases() {
		return new String(getReadBases(), getDirection() == BreakpointDirection.Backward ? 0 : getAnchorLength(), getBreakpointLength(), StandardCharsets.US_ASCII);
	}
	public String getAnchorBases() {
		return new String(getReadBases(), getDirection() == BreakpointDirection.Forward ? 0 : getBreakpointLength(), getAnchorLength(), StandardCharsets.US_ASCII);
	}
	public int getAnchorLength() {
		if (getDirection() == BreakpointDirection.Forward) return getCigar().getCigarElement(0).getLength();
		else return getCigar().getCigarElement(1).getLength(); 
	}
	public int getBreakpointLength() {
		if (getDirection() == BreakpointDirection.Forward) return getCigar().getCigarElement(1).getLength();
		else return getCigar().getCigarElement(0).getLength(); 
	}
	public float getAssemblyQuality() {
		int quality = Integer.MAX_VALUE;
		byte[] qual = getBaseQualities(); 
		for (int i = 0; qual != null && i < qual.length; i++) {
			quality = Math.min(quality, qual[i]);
		}
		if (quality == Integer.MAX_VALUE) return 0;
		return quality;
	}
	/**
	 * Number of reads in assembly
	 * @return
	 */
	public Integer getReadCount() {
		return (Integer)getAttribute(SocratesSamTags.ASSEMBLER_READ_COUNT);
	}
}
