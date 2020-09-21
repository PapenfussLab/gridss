#include "fastq.h"

static int FASTQ_PHRED_QUAL_OFFSET = 33;
static int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

int write_fastq(bam1_t *record, BGZF *out, kstring_t *linebuf, int start, int end, int unique_name) {
	int i;
	// reset buffer to start
	linebuf->l = 0;
	ks_resize(linebuf, 6 + record->core.l_qname + 2 * record->core.l_qseq + (unique_name ? 2 : 0));
	// write read name
	kputc('@', linebuf);
	kputs(bam_get_qname(record), linebuf);
	if (unique_name && (record->core.flag & BAM_FPAIRED)) {
		if (record->core.flag & BAM_FREAD1) {
			kputs("/1", linebuf);
		} else if (record->core.flag & BAM_FREAD2) {
			kputs("/2", linebuf);
		}
	}
	kputc('\n', linebuf);
	if (bam_is_rev(record)) {
		for (i = end - 1; i >= start; i--) {
			kputc(seq_nt16_str[seq_comp_table[bam_seqi(bam_get_seq(record), i)]], linebuf);
		}
	} else {
		for (i = start; i < end; i++) {
			kputc(seq_nt16_str[bam_seqi(bam_get_seq(record), i)], linebuf);
		}
	}
	kputs("\n+\n", linebuf);
	char *qual = bam_get_qual(record);
	if (qual && *qual != '\xff') {
		size_t qual_len = strnlen(qual, record->core.l_qseq);
		if (end <= qual_len) {
			if (bam_is_rev(record)) {
				for (i = end - 1; i >= start; i--) {
					kputc(FASTQ_PHRED_QUAL_OFFSET + qual[i], linebuf);
				}
			} else {
				for (i = start; i < end; i++) {
					kputc(FASTQ_PHRED_QUAL_OFFSET + qual[i], linebuf);
				}
			}
		}
	}
	kputc('\n', linebuf);
	return bgzf_write(out, ks_str(linebuf), ks_len(linebuf));
}