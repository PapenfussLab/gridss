#ifndef FASTQC_H
#define FASTQC_H

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

int write_fastq(bam1_t *record, BGZF *out, kstring_t *linebuf, int start, int end, int unique_name);

#endif