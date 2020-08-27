#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include "fastq.h"

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   gridss_cport  main_extractFragmentsToFastq [options] [input.bam]\n\n");
	fprintf(stderr, "Options: -o FILE output to file\n");
	fprintf(stderr, "		 -m INT   Minimum length of sequence export.Minimum length of sequence export. [20]\n");
	fprintf(stderr, "		 -x	   exclude soft clipped bases.\n");
	fprintf(stderr, "		 -u	   Ensure exported names are unique by suffixing with '/1' or '/2'\n");
	fprintf(stderr, "		 -@	   Number of threads to use\n");
	fprintf(stderr, "\n");
	return EXIT_FAILURE;
}

static int process(
		BGZF *out,
		kstring_t *linebuf,
		int min_sequence_length,
		int include_soft_clips,
		int unique_names,
		bam1_t *record) {
	if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
		// Only process primary alignments
		return 0;
	}
	int start = 0, end = 0;
	if (record->core.flag & BAM_FUNMAP) {
		// unmapped = full length
		start = 0;
		end = record->core.l_qseq - 1;
	} else if (include_soft_clips) {
		// check for start/end soft clip
		int start_clip_length = 0;
		int end_clip_length = 0;
		uint32_t* cigar = bam_get_cigar(record);
		int n_cigar = record->core.n_cigar;
		for (int i = 0; i < n_cigar; i++) {
			if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
				start_clip_length += bam_cigar_oplen(cigar[i]);
			} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
			} else {
				break;
			}
		}
		for (int i = n_cigar - 1; i >= 0; i--) {
			if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
				end_clip_length += bam_cigar_oplen(cigar[i]);
			} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
			} else {
				break;
			}
		}
		// #######################
		// TODO: parse SA tag to determine the longest unmapped sequence
		if (start_clip_length >= min_sequence_length && start_clip_length > end_clip_length) {
			start = 0;
			end = start_clip_length;
		} else if (end_clip_length > min_sequence_length) {
			end = record->core.l_qseq - 1;
			start = end - end_clip_length + 1;
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
	if (argc != optind + 1) {
		return usage();
	}
	// intialisation
	fp_in = hts_open(argv[optind], "r");
	if (!fp_in) {
		fprintf(stderr, "Unable to open input file %s.\n", argv[optind]);
		goto error;
	}
	hdr = sam_hdr_read(fp_in);
	if (hdr == NULL) {
		fprintf(stderr, "Unable to read SAM header.\n");
		goto error;
	}
	if (!hts_set_threads(fp_in, threads)) {
		fprintf(stderr, "Unable to initialise %d threads for reading input", threads);
		goto error;
	}
	fp_out = bgzf_open(out_filename, "wu");
	if (!fp_out) {
		fprintf(stderr, "Unable to open output file");
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
cleanup:
	ks_free(&linebuf);
	if (record) bam_destroy1(record);
	if (hdr) bam_hdr_destroy(hdr);
	if (fp_in) hts_close(fp_in);
	if (fp_out) bgzf_close(fp_out);
	return status;
error:
	status = EXIT_FAILURE;
	perror("main_unmappedSequencesToFastq failed");
	goto cleanup;
}