#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o xtrace


# Will test following variants
# Y       21154600        C       A       1.399104        12.78
# Y       21154600        C       T       1.346947        12.51

dest_file='test/test_cmd.tsv'
test_data='test/exac_cadd.tsv.gz'

# First variant set
./tabix_filter ${test_data} -v 4 'Y:21154600;A' 'Y:21154600;T' > ${dest_file}

cadd_score=`head -1 ${dest_file} | cut -f 6`
test ${cadd_score} = '12.78'

cadd_score=`tail -1 ${dest_file} | cut -f 6`
test ${cadd_score} = '12.51'

rm ${dest_file}
