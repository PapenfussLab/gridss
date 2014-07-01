package au.edu.wehi.idsv.debruijn.windowed;

import htsjdk.samtools.SAMRecord;

import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceEndCoordinateComparator;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadEvidenceAssembler;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.SAMRecordUtil;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.sam.AnomolousReadAssembly;

import com.google.common.collect.Lists;

/**
 * Generates local breakpoint read de bruijn graph assemblies of SV-supporting reads
 * in a single coordinate-sorted pass over the read evidence.
 * 
 * An assembly is generated for each non-reference contig.
 * 
 * @author Daniel Cameron
 *
 */
public class DeBruijnWindowedAssembler implements ReadEvidenceAssembler {
	// TODO: special case breakends that we can assemble across forward and reverse graphs
	// to create a breakpoint
	// Should also genotype the events (INS, DEL, INV, DUP, ...)	
	/**
	 * INS:
	 * contig   ______
	 *          \    /
	 * ref    ___\  /___
	 */
	/**
	 * DEL:
	 *            /\
	 * ref    ___/  \___
	 */
}
