#include <stdio.h>
#include <unistd.h>
#include <htslib/khash.h>
#include "fastq.h"

#define MAX_READNAME_LENGTH 255

static khint_t bam_readname_hash_func(bam1_t* r) {
	const char* str = bam_get_qname(r);
	return kh_str_hash_func(str);
}
static khint_t bam_readname_hash_equal(bam1_t* a, bam1_t* b) {
	const char* stra = bam_get_qname(a);
	const char* strb = bam_get_qname(b);
	return kh_str_hash_equal(stra, strb);
}

KHASH_SET_INIT_STR(str)
KHASH_INIT(bam_readname, bam1_t*, char, 0, bam_readname_hash_func, bam_readname_hash_equal)

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   gridss_cport extractFragmentsToFastq [options] -r <readnames.txt> <input.bam>\n\n");
	fprintf(stderr, "Options: -o FILE  fastq output for unpaired reads.\n");
	fprintf(stderr, "         -1 FILE  fastq output for first in pair reads.\n");
	fprintf(stderr, "         -2 FILE  fastq output for second in pair reads.\n");
	fprintf(stderr, "         -r FILE  File contain names of reads. 1 per line.");
	fprintf(stderr, "         -@       Number of threads to use\n");
	fprintf(stderr, "\n");
	return EXIT_FAILURE;
}

static void output(
		BGZF *out_unpaired,
		BGZF *out_first,
		BGZF *out_second,
		kstring_t *linebuf,
		bam1_t *record1,
		bam1_t *record2) {
	if (record2 && record1) {
		if (out_first) write_fastq(record1, out_first, linebuf, 0, record1->core.l_qseq, 0);
		if (out_second) write_fastq(record2, out_second, linebuf, 0, record2->core.l_qseq, 0);
	} else if (record1) {
		if (out_unpaired) write_fastq(record1, out_unpaired, linebuf, 0, record1->core.l_qseq, 0);
	} else if (record2) {
		if (out_unpaired) write_fastq(record2, out_unpaired, linebuf, 0, record2->core.l_qseq, 0);
	}
}

// returns 0 is this fragment is considered to be processed
static int process(
		BGZF *out_unpaired,
		BGZF *out_first,
		BGZF *out_second,
		kstring_t *linebuf,
		khash_t(bam_readname) *mate_lookup,
		bam1_t *record) {
	if (record->core.flag & BAM_FPAIRED) {
		khiter_t it = kh_get(bam_readname, mate_lookup, record);
		if (it == kh_end(mate_lookup)) {
			int ret;
			bam1_t* inserted_record = bam_dup1(record);
			kh_put(bam_readname, mate_lookup, inserted_record, &ret);
			return 1;
		} else {
			bam1_t *mate = kh_key(mate_lookup, it);
			output(out_unpaired, out_first, out_second, linebuf,
				(record->core.flag & BAM_FREAD1) ? record : mate,
				(record->core.flag & BAM_FREAD1) ? mate : record);
			kh_del(bam_readname, mate_lookup, it);
			bam_destroy1(mate);
			return 0;
		}
	} else {
		output(out_unpaired, out_first, out_second, linebuf, record, NULL);
		return 0;
	}
}

static khash_t(str)* load_readnames(char* filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) return NULL;
	size_t buffersize = MAX_READNAME_LENGTH + 1;
	char * buffer = malloc(sizeof(char) * buffersize);
	khash_t(str) *lookup = kh_init(str);
	khint_t k;
	int absent;
	int line_length = 0;
	while ((line_length = getline(&buffer, &buffersize, fp)) != -1) {
		if (buffer[line_length - 1] == '\n') {
			buffer[line_length - 1] = '\0';
		}
		k = kh_put(str, lookup, strndup(buffer, MAX_READNAME_LENGTH + 1), &absent);
	}
	free(buffer);
	fclose(fp);
	return lookup;
}

int main_extractFragmentsToFastq(int argc, char *argv[]) {
	char *outu_filename = NULL, *out1_filename = NULL, *out2_filename = NULL, *readname_filename = NULL;
	int threads = 1;
	htsFile *fp_in = NULL;
	BGZF *fp_outu = NULL, *fp_out1 = NULL, *fp_out2 = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *record = bam_init1();
	int status = EXIT_SUCCESS;
	kstring_t linebuf = KS_INITIALIZE;
	// arg parsing
	int c;
	while ((c = getopt(argc, argv, "o:1:2:r:@:")) >= 0) {
		switch (c) {
			case 'o': outu_filename = strdup(optarg); break;
			case '1': out1_filename = strdup(optarg); break;
			case '2': out2_filename = strdup(optarg); break;
			case 'r': readname_filename = strdup(optarg); break;
			case '@': threads = strtol(optarg, 0, 0); break;
			default: return usage();
		}
	}
	if (argc != optind + 1) {
		return usage();
	}
	if (!readname_filename) {
		usage();
		fprintf(stderr, "\n-r is a required argument.\n");
		goto error;
	}
	if (!outu_filename && !out1_filename && !out2_filename) {
		usage();
		fprintf(stderr, "\nAt least one of -o -1 -2 must be specified for any output to be produced.\n");
		goto error;
	}
	khash_t(str) *readname_lookup = NULL;
	if (!(readname_lookup = load_readnames(readname_filename))) {
		fprintf(stderr, "Unable to load read name file %s\n", readname_filename);
		goto error;
	}
	fprintf(stderr, "Found %d distinct read name in %s\n", kh_size(readname_lookup), readname_filename);
	char* in_filename = argv[optind];
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
	if (outu_filename) {
		if (!(fp_outu = bgzf_open(outu_filename, "wu"))) {
			fprintf(stderr, "Unable to open output file %s\n", outu_filename);
			goto error;
		}
	}
	if (out1_filename) {
		if (!(fp_out1 = bgzf_open(out1_filename, "wu"))) {
			fprintf(stderr, "Unable to open output file %s\n", out1_filename);
			goto error;
		}
	}
	if (out2_filename) {
		if (!(fp_out2 = bgzf_open(out2_filename, "wu"))) {
			fprintf(stderr, "Unable to open output file %s\n", out2_filename);
			goto error;
		}
	}
	// process records
	khash_t(bam_readname) *mate_lookup = kh_init(bam_readname);
	int r;
	khiter_t it;
	while ((r = sam_read1(fp_in, hdr, record)) >= 0) {
		if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
		it = kh_get(str, readname_lookup, bam_get_qname(record));
		if (it != kh_end(readname_lookup)) {
			if (!process(fp_outu, fp_out1, fp_out2, &linebuf, mate_lookup, record)) {
				kh_del(str, readname_lookup, it);
			}
		}
	}
	if (r < -1) {
		fprintf(stderr, "Truncated file.\n");
		goto error;
	}
	if (kh_size(mate_lookup) > 0) {
		fprintf(stderr, "Found %d paired records without mates. Treating as unpaired.\n", kh_size(mate_lookup));
		for (khint_t it = 0; it < kh_end(mate_lookup); it++) {
			if (kh_exist(mate_lookup, it)) {
				bam1_t *orphan = kh_key(mate_lookup, it);
				output(fp_outu, fp_out1, fp_out2, &linebuf, orphan, NULL);
				khint_t rnit = kh_get(str, readname_lookup, bam_get_qname(orphan));
				if (rnit != kh_end(readname_lookup)) {
					kh_del(str, readname_lookup, rnit);
				}
				bam_destroy1(orphan);
			}
		}
		kh_destroy(bam_readname, mate_lookup);
	}
	if (kh_size(readname_lookup) > 0) {
		fprintf(stderr, "Could not find primay alignment records for %d read names.\n", kh_size(readname_lookup));
	}
cleanup:
	if (fp_outu) bgzf_close(fp_outu);
	if (fp_out1) bgzf_close(fp_out1);
	if (fp_out2) bgzf_close(fp_out2);
	if (readname_lookup) {
		for (khint_t it = 0; it < kh_end(readname_lookup); it++) {
			if (kh_exist(readname_lookup, it)) {
				free((char*)kh_key(readname_lookup, it));
			}
		}
		kh_destroy(str, readname_lookup);
	}
	ks_free(&linebuf);
	if (record) bam_destroy1(record);
	if (hdr) bam_hdr_destroy(hdr);
	if (fp_in) hts_close(fp_in);
	return status;
error:
	status = EXIT_FAILURE;
	perror("main_extractFragmentsToFastq failed");
	goto cleanup;
}