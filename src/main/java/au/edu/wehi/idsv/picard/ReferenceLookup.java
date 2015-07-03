package au.edu.wehi.idsv.picard;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public interface ReferenceLookup extends ReferenceSequenceFile {
	public byte getBase(int referenceIndex, int position);
}
