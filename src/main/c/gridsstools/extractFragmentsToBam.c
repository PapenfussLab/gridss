#include <stdio.h>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/thread_pool.h>

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
	fprintf(stderr, "Usage:   gridsstools extractFragmentsToBam [options] <readnames.txt> <input.bam>\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "         -@       Number of threads to use\n");
	fprintf(stderr, "         -o       Output file\n");
	fprintf(stderr, "\n");
	return EXIT_FAILURE;
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

int main_extractFragmentsToBam(int argc, char *argv[]) {
	int threads = 1;
	htsFile *fp_in = NULL;
	htsFile *fp_out = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *record = bam_init1();
	int status = EXIT_SUCCESS;
	char *strargs = stringify_argv(argc, argv);
	char *out_filename = "-";
	// arg parsing
	int c;
	while ((c = getopt(argc, argv, "o:@:")) >= 0) {
		switch (c) {
			case 'o': out_filename = strdup(optarg); break;
			case '@': threads = strtol(optarg, 0, 0); break;
			default: return usage();
		}
	}
	if (argc != optind + 2) {
		return usage();
	}
	char *readname_filename = argv[optind];
	char *in_filename = argv[optind + 1];
	khash_t(str) *readname_lookup = NULL;
	if (!(readname_lookup = load_readnames(readname_filename))) {
		fprintf(stderr, "Unable to load read name file %s\n", readname_filename);
		goto error;
	}
	fprintf(stderr, "Read %d distinct read names\n", kh_size(readname_lookup), readname_filename);
	if (!(fp_in = hts_open(in_filename, "r"))) {
		fprintf(stderr, "Unable to open input file %s.\n", in_filename);
		goto error;
	}
	htsThreadPool p = {NULL, 0};
	if (threads > 0) {
		if (!(p.pool = hts_tpool_init(threads))) {
			fprintf(stderr, "Error creating thread pool\n");
			goto error;
		}
		hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &p);
	}
	hdr = sam_hdr_read(fp_in);
	if (hdr == NULL) {
		fprintf(stderr, "Unable to read SAM header.\n");
		goto error;
	}
	if (!(fp_out = hts_open(out_filename, "w"))) {
		fprintf(stderr, "Unable to open output file %s.\n", out_filename);
		goto error;
	}
	hts_set_opt(fp_out, HTS_OPT_THREAD_POOL, &p);
	if (sam_hdr_write(fp_out, hdr) != 0) {
		fprintf(stderr, "Unable to write SAM header.\n");
		goto error;
	}
	// process records
	int r;
	khiter_t it;
	while ((r = sam_read1(fp_in, hdr, record)) >= 0) {
		// could switch to a hash-based lookup for long read names
		it = kh_get(str, readname_lookup, bam_get_qname(record));
		if (it != kh_end(readname_lookup)) {
			if (sam_write1(fp_out, hdr, record) < 0) {
				fprintf(stderr, "Error writing %s to output file.\n", bam_get_qname(record));
				goto error;
			}
		}
	}
	if (r < -1) {
		fprintf(stderr, "Truncated file.\n");
		goto error;
	}
cleanup:
	if (fp_in) hts_close(fp_in);
	if (fp_out) hts_close(fp_out);
	if(strargs) free(strargs);
	if (record) bam_destroy1(record);
	if (hdr) bam_hdr_destroy(hdr);
	if (readname_lookup) {
		for (it = 0; it < kh_end(readname_lookup); it++) {
			if (kh_exist(readname_lookup, it)) {
				free((char*)kh_key(readname_lookup, it));
			}
		}
		kh_destroy(str, readname_lookup);
	}
	return status;
error:
	status = EXIT_FAILURE;
	perror("main_extractFragmentsToFastq failed");
	goto cleanup;
}