package au.edu.wehi.socrates.debruijn;

import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import au.edu.wehi.socrates.BreakpointDirection;
import au.edu.wehi.socrates.BreakpointLocation;
import au.edu.wehi.socrates.DirectedBreakpointAssembly;
import au.edu.wehi.socrates.DirectedEvidence;
import au.edu.wehi.socrates.DirectedEvidenceEndCoordinateComparator;
import au.edu.wehi.socrates.LinearGenomicCoordinate;
import au.edu.wehi.socrates.NonReferenceReadPair;
import au.edu.wehi.socrates.ReadEvidenceAssembler;
import au.edu.wehi.socrates.SoftClipEvidence;
import au.edu.wehi.socrates.sam.AnomolousReadAssembly;

import com.google.common.collect.Lists;
import com.google.common.collect.Queues;


