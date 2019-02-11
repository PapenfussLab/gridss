#!/usr/bin/env python

import os
import sys

BLACKLIST_FILENAME = "data/ENCODE_blacklist_hg19/ENCFF001TDO.bed"
REFERENCE_GENOME = "data/hg19/hg19.fa"
GRIDSS_JARFILE = "./gridss-2.2.0-gridss-jar-with-dependencies.jar"


# Read the script template
with open("somatic_template.sh") as template_file:
    script_template = template_file.read()


# Take the patient sampledata filename from the command line
sample_file_list = sys.argv[1]

# Parse the patient sample data
with open(sample_file_list) as infile, open("3_run_scripts.sh", "w") as runfile:
    for line in infile:
        line = line.strip()
        if line[0]=="#": continue
        print(line)
        tokens = line.split(", ")

        patient_id = tokens[0]

        inputs_string = []
        for token in tokens[3:]:
            inputs_string.append('INPUT="%s"' % token)
        inputs_string = " ".join(inputs_string)

        parameters_dict = {
            "patient_id": patient_id,
            "somatic_sv_vcf": tokens[1],
            "output_assembly_bam": tokens[2],
            "inputs": inputs_string,
            "blacklist_file": BLACKLIST_FILENAME,
            "reference_genome": REFERENCE_GENOME,
            "gridss_jarfile": GRIDSS_JARFILE
        }

        script = script_template % parameters_dict

        outfilename = "scripts/%s.sh" % patient_id
        with open(outfilename, "w") as outfile:
            print(script, file=outfile)
            os.system("chmod 755 %s" % outfilename)

        runfile.write(outfilename + "\n")
