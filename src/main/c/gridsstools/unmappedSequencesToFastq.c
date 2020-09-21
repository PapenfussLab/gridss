#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <ctype.h>
#include "fastq.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   gridss_cport  extractFragmentsToFastq [options] [input.bam] ...\n\n");
	fprintf(stderr, "Options: -o FILE  output to file\n");
	fprintf(stderr, "         -m INT   Minimum length of sequence export.Minimum length of sequence export. [20]\n");
	fprintf(stderr, "         -x       exclude soft clipped bases.\n");
	fprintf(stderr, "         -u       Ensure exported names are unique by suffixing with '/1' or '/2'\n");
	fprintf(stderr, "         -@       Number of threads to use\n");
	fprintf(stderr, "\n");
	return EXIT_FAILURE;
}

static int parse_cigar_op(
		char** str,
		int* op,
		int* len) {
	*len = 0;
	while (isdigit(**str)) {
		*len *= 10;
		*len += **str - '0';
		(*str)++;
	}
	if (len == 0) {
		return 0;
	}
	*op = bam_cigar_table[(unsigned char)**str];
	(*str)++;
	return *op != -1;
}

// Extracts the start and end clipping lengths from the chimeric alignment in the given SA tag
// str: [in/out] parameter of current position in the SA string.
// 		output value is the position after the SA chimeric alignment,
// 		or the position the parsing failure
// start_clip_length: [out] length of start clipping (hard and soft)
// end_clip_length: [out] length of end clipping (hard and soft)
// returns true if the parsing was successful
static int parse_SA_tag_clipping(
		char** str,
		int* start_clip_length,
		int* end_clip_length,
		int* is_negative_strand) {
	// SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
	*start_clip_length = 0;
	*end_clip_length = 0;
	*is_negative_strand = 0;
	*str = strchr(*str, ','); if (!*str) return 0;
	*str += 1;
	*str = strchr(*str, ','); if (!*str) return 0;
	*str += 1;
	*is_negative_strand = **str == '-';
	*str = strchr(*str, ','); if (!*str) return 0;
	*str += 1;
	// process CIGAR
	int oplen;
	int op;
	int inStart = 1;
	while (parse_cigar_op(str, &op, &oplen)) {
		if (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) {
			if (inStart) {
				*start_clip_length += oplen;
			} else {
				*end_clip_length += oplen;
			}
		} else {
			inStart = 0;
		}
	}
	// advance to end of chimeric alignment record
	while (*str && **str != ';' && **str != '\0') {
		(*str)++;
	}
	if (*str && **str == ';') (*str)++;
	return 1;
}

static int process(
		BGZF *out,
		kstring_t *linebuf,
		int min_sequence_length,
		int include_soft_clips,
		int unique_names,
		bam1_t *record) {
	int i;
	if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
		// Only process primary alignments
		return 0;
	}
	int start = 0, end = 0;
	if (record->core.flag & BAM_FUNMAP) {
		// unmapped = full length
		start = 0;
		end = record->core.l_qseq;
	} else if (include_soft_clips) {
		// check for start/end soft clip
		int start_clip_length = 0;
		int end_clip_length = 0;
		uint32_t* cigar = bam_get_cigar(record);
		int n_cigar = record->core.n_cigar;
		for (i = 0; i < n_cigar; i++) {
			if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
				start_clip_length += bam_cigar_oplen(cigar[i]);
			} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
			} else {
				break;
			}
		}
		for (i = n_cigar - 1; i >= 0; i--) {
			if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
				end_clip_length += bam_cigar_oplen(cigar[i]);
			} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
			} else {
				break;
			}
		}
		// SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
		uint8_t *satag = bam_aux_get(record, "SA");
		if (satag) {
			char *sa = bam_aux2Z(satag);
			if (sa) {
				int sa_start_clip_length = 0;
				int sa_end_clip_length = 0;
				int saNegativeStrand;
				while(parse_SA_tag_clipping(&sa, &sa_start_clip_length, &sa_end_clip_length, &saNegativeStrand)) {
					if ((saNegativeStrand && !bam_is_rev(record)) ||
							(!saNegativeStrand && bam_is_rev(record))) {
						// strand differs from our record - need to flip start and end
						int tmp = sa_start_clip_length;
						sa_start_clip_length = sa_end_clip_length;
						sa_end_clip_length = tmp;
					}
					start_clip_length = sa_start_clip_length < start_clip_length ? sa_start_clip_length : start_clip_length;
					end_clip_length = sa_end_clip_length < end_clip_length ? sa_end_clip_length : end_clip_length;
				}
			}
		}
		if (start_clip_length >= min_sequence_length && start_clip_length > end_clip_length) {
			start = 0;
			end = start_clip_length;
		} else if (end_clip_length >= min_sequence_length) {
			end = record->core.l_qseq;
			start = record->core.l_qseq - end_clip_length;
		}
	}
	if (end - start >= min_sequence_length) {
		write_fastq(record, out, linebuf, start, end, unique_names);
	}
}

int main_unmappedSequencesToFastq(int argc, char *argv[]) {
	char *out_filename = "-";
	int min_sequence_length = 20, include_soft_clips = 1, unique_names = 0, threads = 1;
	htsFile *fp_in = NULL;
	BGZF *fp_out = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *record = bam_init1();
	int status = EXIT_SUCCESS;
	kstring_t linebuf = KS_INITIALIZE;
	int i;
	// arg parsing
	int c;
	while ((c = getopt(argc, argv, "o:m:xu@:")) >= 0) {
		switch (c) {
			case 'o': out_filename = strdup(optarg); break;
			case 'm': min_sequence_length = strtol(optarg, 0, 0); break;
			case 'x': include_soft_clips = 0; break;
			case 'u': unique_names = 1; break;
			case '@': threads = strtol(optarg, 0, 0); break;
			default: return usage();
		}
	}
	if (argc < optind + 1) {
		return usage();
	}
	for (i = optind; i < argc; i++) {
		char* in_filename = argv[i];
		// open files
		if (!(fp_in = hts_open(in_filename, "r"))) {
			fprintf(stderr, "Unable to open input file %s.\n", in_filename);
			goto error;
		}
		if (hts_set_threads(fp_in, threads)) {
			fprintf(stderr, "Unable to create thread pool with %d threads\n", threads);
			goto error;
		}
		hdr = sam_hdr_read(fp_in);
		if (hdr == NULL) {
			fprintf(stderr, "Unable to read SAM header.\n");
			goto error;
		}
		fp_out = bgzf_open(out_filename, "wu");
		if (!fp_out) {
			fprintf(stderr, "Unable to open output file\n");
			goto error;
		}
		// process records
		int r;
		while ((r = sam_read1(fp_in, hdr, record)) >= 0) {
			process(fp_out, &linebuf, min_sequence_length, include_soft_clips, unique_names, record);
		}
		if (r < -1) {
			fprintf(stderr, "Truncated file.\n");
			goto error;
		}
		if (hdr) {
			bam_hdr_destroy(hdr);
			hdr = NULL;
		}
		if (fp_in) {
			hts_close(fp_in);
			fp_in = NULL;
		}
	}
cleanup:
	ks_free(&linebuf);
	if (record) bam_destroy1(record);
	if (fp_out) bgzf_close(fp_out);
	if (hdr) bam_hdr_destroy(hdr);
	if (fp_in) hts_close(fp_in);
	return status;
error:
	status = EXIT_FAILURE;
	perror("main_unmappedSequencesToFastq failed");
	goto cleanup;
}