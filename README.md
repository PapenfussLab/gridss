README.txt
Authors: Arthur Hsu, Jan Schroeder, Anthony T Papenfuss 
Date: 27/3/2013

A driver script (developed in Python2.7) - Socrates, can be found in this directory.
For a test example run with the script, download the test data from http://bioinf.wehi.edu.au/socrates/test_data.tar.gz, extract in the script directory and execute:
./Socrates all test_data/data/random_ecoli_63_s.bam test_data/data/bowtie2_db/ecoli1665_bt2

Additional documentation is currently being developed. 

The Socrates package contains several Java programs (developed in JDK1.6).
Each program is designed to process data for a specific stages of analysis.
Below is a list of the included programs and each of them will be elaborated
latter.

1) PreprocessBAM
2) ProcessRealignmentBAM
3) PredictRearrangements
4) AnnotateRearrangements

These are simple wrapper scripts that sets library and execution path to run
Socrates programs under Linux and Mac.
While default values have worked satisfactorily in our simulated and real
cancer genome sequencing datasets, users should set program parameters in the
driver scripts appropri- ately for their own data. Full lists of program
parameters are provided in the following sections, together with a discussion
on the impact of changing them where applicable.


1.1. Preprocess BAM File 
Stratifies the original BAM file into...

usage: PreprocessBAM [options] alignment bam
-b,	--base-quality <score>	Minimum average base quality score of soft clipped sequence  [default: 5] 
-h,	--help 			print this message 
-k,	--keep-duplicate	keep duplicate reads [default: false]
-l,	--long-sc-len <length>	Length threshold of long soft-clip [default: 25 (bp)]
-p,	--percent-id <pid> 	Minimum alignment percent identity to reference [default: 95 (%)]
-q,	--min-mapq <mapq> 	Minimum alignments mapq [default: 5]
-t,	--threads <threads> 	Number of threads to use [default: 1]
-v,	--verbose		be verbose of progress
    
Minimum base quality option: 
A reasonable threshold helps removing low quality
soft clips that could lead to erroneous breakpoint calls.

Long soft clip length: 
Studies have shown that longer the sequences the more
likely they can be uniquely placed in a genome. In an early study, it is
demonstrated that while percentage of unique mapping improves with increasing
read length, the rate of gain di- minishes past 25nt ( 80% at 25nt and 90% at
40nt). If value for this parameter is too low, many non-unique soft clips will
be produced and impact on system requirement, processing time and reliability
of results downstream. On the other hand, too high the value results in low
number of long soft clips and hence risk of missing breakpoints.
Percent identity: We often observe higher-than-expected base mismatch rate for
reads in satellite, centromeric and telomeric regions where correctness of
alignments can be con- tentious. Minimum percent identity threshold, which is
equivalent to maximum allowable mismatch rate, can greatly reduce these
erroneous alignments.

Minimum mapping quality: 
Higher mapping quality, while may not guarantee
unique align- ment, is sufficient to exclude multi-mapping anchor alignments
from further analysis for Bowtie2 and BWA aligned reads.


1.2. Process the re-alignment BAM file
usage: ProcessRealignmentBAM input_bam output_bam anchor_info
input_bam 	Re-aligned soft clip BAM file. Use “-” to accept input from stdin
output_bam 	Output re-alignment BAM with anchor info merged
anchor_info	Anchor info file produced by BAMStratifier

 This program merges soft clip re-alignment BAM file with anchor alignment
information. The program has built-in sorting mechanism and therefore can take
unsorted, raw re- alignment output from aligner. While the program accepts
input BAM file from standard input channel, this requires more system memory
for buffering.


1.3. Predict rearrangements
usage: PredictRearrangements [options] realigned_sc_bam short_sc_bam metrics_file
-f,	--flank <flank> 		Size of flank for promiscuity filter [default: 50 (bp)]
-h,	--help 				print this message 
-i,	--ideal-only 			Use only proper pair 5’ SC and anomalous pair 3’ SC [default: false]
-l,	--long-sc-len <length> 		Length threshold of long soft-clip [default: 25 (bp)]
-m,	--promiscuity <threshold>	Exclude cluster if more than [promiscuity] clusters within [flank]bp of a break
-p,	--percent-id <pid>		Minimum realignment percent identity to reference [default: 95 (%)]
-q,	--min-mapq <mapq> 		Minimum realignments mapq [default: 2]
-c,	--short-sc-cluster 		Search for short soft clip cluster support for unpaired clusters 
-s,	--max-support <support> 	Maximum realignment support to search for short SC cluster [default: 30]
-t,	--threads <threads> 		Number of threads to use [default: 3]
-v,	--verbose			be verbose of progress
     
  
1.4. Annotating rearrangements
usage: AnnotateRearrangements [options] socrates_paired_cluster_output
-n,	--normal <normal>	Socrates paired breakpoint calls for normal sample 
-r,	--repeatmask <file>	UCSC repeat masker track file in BED format, Tabix indexed.


