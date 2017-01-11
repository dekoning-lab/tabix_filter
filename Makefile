CFLAGS=-Wall -std=c99
ASANFLAGS=-g -O1 -fsanitize=address -fno-omit-frame-pointer
INCPATH=./include
LIBPATH=/usr/local/lib
LIBS=-lhts -lm

all: tabix_filter

tabix_filter: tabix_filter.c
		${CC} $^ ${CFLAGS} -I${INCPATH} -L${LIBPATH} ${LIBS} -o $@

tabix_filter_debug_osx: tabix_filter.c
		${CC} $^ ${CFLAGS} ${ASANFLAGS} -I${INCPATH} -L${LIBPATH} ${LIBS} -o $@

clean:
		rm -rf tabix_filter tabix_filter.dSYM

CADD_URL='http://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3.tsv.gz'

test: tabix_filter test/exac_cadd.tsv.gz test/exac_cadd.tsv.gz.tbi test_cmd test_file

test/exac_cadd.tsv.gz test/exac_cadd.tsv.gz.tbi:
	curl ${CADD_URL} -o test/exac_cadd.tsv.gz
	tabix -s1 -b2 -e2 test/exac_cadd.tsv.gz

test_cmd: test/exac_cadd.tsv.gz
	sh test/test_cmd.sh

test_file: test/exac_cadd.tsv.gz
	sh test/test_file.sh

test_clean:
	rm test/exac_cadd.tsv.gz test/exac_cadd.tsv.gz

