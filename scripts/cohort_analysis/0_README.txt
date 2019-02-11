README.txt
==========

Author: Tony Papenfuss
Email: papenfuss@wehi.edu.au


Summary
-------

GRIDSS analyses multiple bams from a single patient simulataneously, thus
achieving high sensitivity and consistency between samples. This software
package provides helper scripts to easily run GRIDSS on bam files from a
cohort of patients with one or more samples per patient.

To run GRIDSS on a cohort with potentially different numbers of samples per
patient. Each line provides the patient id, output VCF and bam files, and
a variable number of input bam files (germline, primary, regional, distant met,
...).

This readme guides the user through the installation and use of GRIDSS for this.


Dependencies: java 1.8, bwa mem

GRIDSS source code and documentation is available from https://github.com/PapenfussLab/gridss.
An example script to run GRIDSS is available in https://github.com/PapenfussLab/gridss/blob/master/example/somatic.sh


Instructions
------------

1. Check that java 1.8 and bwa mem are installed.


2. On the command line run

> ./1_setup.sh

This will download the GRIDSS jar file and dependencies
(https://github.com/PapenfussLab/gridss/releases/download/v2.2.0/gridss-2.2.0-gridss-jar-with-dependencies.jar)
and download and decompress the ENCODE blacklist file
(https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz).


3. If necessary index your reference genome:

> samtool faidx <reference_genome_location>
> bwa index <reference_genome_location>


4. Open 2_generate_analysis_scripts.py using an editor and edit the locations
of the following files:

BLACKLIST_FILENAME = "data/ENCODE_blacklist_hg19/ENCFF001TDO.bed"
REFERENCE_GENOME = "data/hg19.fa"
GRIDSS_JARFILE = "./gridss-2.2.0-gridss-jar-with-dependencies.jar"


5. Open & edit the samples file sample.csv using a text editor. The file
contains one line per patient. All related samples are processed by GRIDSS
together. The columns of this file are as follows:

- patient_id
- GRIDSS output somatic SV vcf file
- GRIDSS output assembly bam
- normal bam
- first tumour bam (typically primary)
- as many additional comma deparated bams as needed


6. Run

> ./2_generate_analysis_scripts.py sample.csv

This will generate scripts for running GRIDSS and a file (3_run_scripts.sh) for
running these.


7. Run

> ./3_run_scripts.sh

or the scripts individually.

8. If you want, cleanup by running:

> ./4_cleanup.sh
