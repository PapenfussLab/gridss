# GRIDSS Docker


## Summary

GRIDSS analyses multiple bams from a single patient simulataneously, thus
achieving high sensitivity and consistency between samples. This software
package provides docker images and helper scripts to easily run GRIDSS on bam files from a
cohort of patients with one or more samples per patient.

To run GRIDSS on a cohort with potentially different numbers of samples per
patient. Each line provides the patient id, output VCF and bam files, and
a variable number of input bam files (germline, primary, regional, distant met,
...).

This readme guides the user through the installation and use of GRIDSS for this.

## Installation

The docker image can be pulled 

## Running

### Create cohort CSV

Open & edit the samples file sample.csv using a text editor. The file
contains one line per patient. All related samples are processed by GRIDSS
together. The columns of this file are as follows:

- patient_id
- GRIDSS output somatic SV vcf file. This can be left blank and a name based on the patient_id will be used
- GRIDSS output assembly bam. This can be left blank and a name based on the patient_id will be used
- normal bam
- first tumour bam (typically primary)
- as many additional comma deparated bams as needed

### Set reference genome location

- Edit `run_cohort_from_csv.sh` and update the `REFERENCE` genome location.

If the default ENCODE blacklist or a genome other than hg19 is used, the `BLACKLIST` location should be updated. An empty BED file can be used if no blacklist is desired.

### Generate per-sample docker scripts

Execute `run_cohort_from_csv.sh sample.csv output_dir` to generate a GRIDSS run script for each sample in `output_dir`.

Note that if you do not have a bwa index for your reference genome, this script will attempt to run bwa from within the GRIDSS docker container to generate one.

### Execute run scripts

Once the run scripts have been generated, these can be executed directly or scheduled on a cluster. In dockerised form, GRIDSS requires up to 24Gb (16 for GRIDSS, 8 for bwa) of memory and 4 cores.