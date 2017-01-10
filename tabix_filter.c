#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <getopt.h>
#include "include/tbx.h"
#include "include/kstring.h"

#define KS_INIT {0, 0, NULL}
#define BUFSIZE 2048
#define KS_ALLOC (kstring_t *)(calloc(1, sizeof(kstring_t)))

#define FILE_NOT_FOUND 2


typedef struct {
    kstring_t *input_path;
    kstring_t *variant_path;
    int var_column;
} args_t;

typedef struct {
    kstring_t *chromosome;
    unsigned long long position;
    kstring_t *variant;
} query_t;

static struct option long_options[] = {
    {"input", required_argument, NULL, 'i'},
    {"variants", required_argument, NULL, 'V'},
    {"var-column", required_argument, NULL, 'H'},
    {NULL, 0, NULL, 0}
};

void usage() {
    fprintf(stderr, "\n"
            "Usage: tabix_variant [OPTIONS] <FILE> [QUERY [...]]\n"
            "\nQuery format:\n"
            "  'CHROM:POS;VAR' e.g. 'Y:2655034;C'\n"
            "\nOptions:\n"
            "  -V --variants FILE       Restrict to variants listed in the file\n"
            "                           Variants file should be a whitespace-delimeted file with CHROM, POS, VAR columns\n"
            "                           Alternatively, location of VAR column can be set with `--var-column` argument\n"
            "  -v --var-column INT      1-indexed location of column with the variant to be retrieved\n");
}

void error(const char *format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}

kstring_t *make_tabix_region(kstring_t *chromosome, unsigned long long position) {
    kstring_t *region = KS_ALLOC;
    ksprintf(region, "%s:%llu-%llu", chromosome->s, position, position);
    return region;
}

size_t count_lines(char *file_name) {
  FILE *f = fopen(file_name, "r");
  if (!f) {
      error("Could not open variant file %s\n", file_name);
      exit(FILE_NOT_FOUND);
  }
  char ch;
  size_t n_lines = 0;
  uint64_t pc = EOF;
  while ((ch = fgetc(f)) != EOF) {
    if (ch == '\n') n_lines++;
    pc = ch;
  }
  if (pc != (uint64_t)EOF && pc != (uint64_t)('\n')) n_lines++;
  fclose(f);
  return n_lines;
}

query_t *parse_queries_cmd(size_t argc, char *argv[argc]) {
    query_t *queries = (query_t *)calloc(argc, sizeof(query_t));
    size_t i = 0;
    for (; i < argc; i++) {

        ks_tokaux_t aux;
        char *p = kstrtok(argv[i], ":", &aux);
        if (p == NULL) goto error;
        queries[i].chromosome = KS_ALLOC;
        kputsn(p, aux.p - p, queries[i].chromosome);

        p = kstrtok(NULL, ";", &aux);
        if (p == NULL) goto error;
        kstring_t *tmp = KS_ALLOC;
        kputsn(p, aux.p - p, tmp);
        queries[i].position = strtoll(tmp->s, NULL, 10);
        free(tmp);

        p = kstrtok(NULL, NULL, &aux);
        if (p == NULL) goto error;
        queries[i].variant = KS_ALLOC;
        kputsn(p, aux.p - p, queries[i].variant);
    }
    // for (int i = 0; i < argc; i++) {
    //     printf("%s:%llu;%s\n", queries[i].chromosome->s, queries[i].position, queries[i].variant->s);
    // }
    return queries;
error:
    error("Can not parse query %s\n", argv[i]);
    free(queries);
    usage();
    exit(EXIT_FAILURE);
}

query_t *parse_queries_file(size_t n_queries, kstring_t *variant_path) {
    query_t *queries = (query_t *)calloc(n_queries, sizeof(query_t));
    if (strcasecmp(variant_path->s + variant_path->l - 4, ".tsv") != 0) {
        error("The variant file %s does not appear to be a .tsv file\n", variant_path->s);
        exit(EXIT_FAILURE);
    }
    FILE *variant_file = fopen(variant_path->s, "r");
    if (!variant_file) {
        error("Could not open variant file %s\n", variant_path);
        exit(FILE_NOT_FOUND);
    }

    int n_columns;
    kstring_t *line_buffer = KS_ALLOC;
    for(size_t i = 0; i < n_queries; i++) {
        line_buffer->l = 0;

        int read = kgetline(line_buffer, (kgets_func *)fgets, variant_file);
        if (read == EOF) break;

        int n_splits;
        int *offsets = ksplit(line_buffer, 0, &n_splits);

        if (i == 0) n_columns = n_splits;
        if (n_columns != n_splits) {
            error("Expected %d columns in line %d, got %d\n", n_columns, i+1, n_splits);
            free(queries);
            exit(EXIT_FAILURE);
        }

        queries[i].chromosome = KS_ALLOC;
        kputs(line_buffer->s + offsets[0], queries[i].chromosome);

        kstring_t *tmp = KS_ALLOC;
        kputs(line_buffer->s + offsets[1], tmp);
        queries[i].position = strtoll(tmp->s, NULL, 10);
        free(tmp);

        queries[i].variant = KS_ALLOC;
        kputs(line_buffer->s + offsets[2], queries[i].variant);

    }
    // for (int i = 0; i < n_queries; i++) {
    //     printf("%s:%llu;%s\n", queries[i].chromosome->s, queries[i].position, queries[i].variant->s);
    // }

    return queries;
}

int get_tabix_lines(char *index_path, size_t len, query_t *queries, size_t column) {

    htsFile *index_file = hts_open(index_path, "r");
    if (!index_file) { 
        error("Could not read %s\n", index_path);
        exit(FILE_NOT_FOUND);
    }

    tbx_t *tabix_index = tbx_index_load(index_path);
    if (!tabix_index) {
        error("Could not load .tbi index of %s\n", index_path);
        exit(FILE_NOT_FOUND);
    }

    kstring_t *str = KS_ALLOC;

    for (size_t i = 0; i < len; i++) {
        kstring_t *region = make_tabix_region(queries[i].chromosome, queries[i].position);
        hts_itr_t *iter = tbx_itr_querys(tabix_index, region->s);
        if (!iter) continue;

        // printf("%s:%llu;%s\n", queries[i].chromosome->s, queries[i].position, queries[i].variant->s);
        while(tbx_itr_next(index_file, tabix_index, iter, str) >= 0) {

            kstring_t *columns = KS_ALLOC;
            kputs(str->s, columns);
            int n_splits;
            int *offsets = ksplit(columns, 0, &n_splits);
            char *variant = columns->s + offsets[column - 1];
            if (strcasecmp(variant, queries[i].variant->s) == 0) {
                printf("%s\n", str->s);
            }
            free(columns);
        }
        tbx_itr_destroy(iter);
        free(region);
    }

    free(str);
    tbx_destroy(tabix_index);
    hts_close(index_file);
    return 0;
}

int main(int argc, char *argv[]) {
    args_t args = {
        .input_path = KS_ALLOC,
        .variant_path = KS_ALLOC,
        .var_column = 3
    };

    while(true) {
        int opt = getopt_long(argc, argv, "i:V:v:", long_options, NULL);
        if (opt == -1) break;

        switch (opt) {
        case 'V':
            kputs(optarg, args.variant_path);
            break;
        case 'v':
            args.var_column = strtol(optarg, NULL, 10);
            break;
        default:
            usage();
            return EXIT_FAILURE;
        }
    }

    query_t *queries = NULL;
    int n_queries = 0;
    if (optind == argc) {
        usage();
        return EXIT_FAILURE;
    } else if (args.variant_path->s == NULL && argc <= optind + 1) {
        error("Must either provide [QUERY] or specify ``--variant_path`\n");
        usage();
        return EXIT_FAILURE;
    } else if (args.variant_path->s != NULL) {
        n_queries = count_lines(args.variant_path->s);
        // printf("parsing file %s with %d lines\n", args.variant_path->s, n_queries);
        queries = parse_queries_file(n_queries, args.variant_path);
    } else if (argc > optind + 1) {
        n_queries = argc - optind - 1;
        queries = parse_queries_cmd(n_queries, argv + optind + 1);
    } else {
        usage();
        return EXIT_FAILURE;
    }

    kputs(argv[optind], args.input_path);
    size_t column = 4;
    get_tabix_lines(args.input_path->s, n_queries, queries, column);

    for (int i = 0; i < n_queries; i++) {
        free(queries[i].chromosome);
        free(queries[i].variant);
    }
    free(queries);
    return 0;
}
