
#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o xtrace

# Will test following variants
#
# Y       21154600        C       A       1.399104        12.78
# Y       21154600        C       T       1.346947        12.51

# Y       25375708        C       A       0.210844        4.793
# Y       25375708        C       T       0.120277        3.825

query_file='test/test_query.tsv'
dest_file='test/test_file.tsv'
test_data='test/exac_cadd.tsv.gz'

cat <<EOF > ${query_file}
Y 21154600 A
Y 21154600 T
Y 25375708 A
Y 25375708 T
EOF

./tabix_filter ${test_data} -v 4 -V ${query_file} > ${dest_file}

cadd_score=`awk 'NR == 1' ${dest_file} | cut -f 6`
test ${cadd_score} = '12.78'

cadd_score=`awk 'NR == 2' ${dest_file} | cut -f 6`
test ${cadd_score} = '12.51'

cadd_score=`awk 'NR == 3' ${dest_file} | cut -f 6`
test ${cadd_score} = '4.793'

cadd_score=`awk 'NR == 4' ${dest_file} | cut -f 6`
test ${cadd_score} = '3.825'

rm ${query_file}
rm ${dest_file}

