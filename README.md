# Tabix-filter

This is a modified version of `tabix` query, which allows to filter SNP by position _and_ by alternative variants.
It can be used as a part of a SNP annotation pipeline which includes multiple variants in the same locus.

## Usage

You can query for specific variant from a command line:

```
tabix annotation.tsv
tabix_filter annotation.tsv 'Y:2655034;C'
```

You can also use a query file for large batch lookups:

```
#variants.tsv
Y 21154600 A
Y 21154600 T
Y 25375708 A
Y 25375708 T
```

```
tabix_filter annotation.tsv -V variants.tsv
```

## Installation

### Dependencies

The code depends on `htslib`.
Please see installation details [here](https://github.com/samtools/htslib).

It should also be possible to install `htslib` with [linuxbrew](http://linuxbrew.sh/) and [homebrew](http://brew.sh/).

```
brew install htslib
```

### Building

Simply run `make`

```bash
make
```

There is also a debug target of MacOS:

```bash
make tabix_filter_debug_osx
```

### Testing

There is a simple set of scripts that test main functionality.

```bash
make test
```

## Full example

This is the example used in `test/test_file.sh`:

```bash
# Download some annotations for CADD:
curl 'http://krishna.gs.washington.edu/download/CADD/v1.3/ExAC_r0.3.tsv.gz' -o exac_cadd.tsv.gz

# Create an index with regular tabix:
tabix -s1 -b2 -e2 exac_cadd.tsv.gz

# Write some data into test_query.tsv file:
cat <<EOF > test_query.tsv
Y 21154600 A
Y 21154600 T
Y 25375708 A
Y 25375708 T
EOF

# Use tabix_filter to retrieve CADD annotations for the variants
# -v specifies which column contains alternative variant
./tabix_filter exac_cadd.tsv.gz -v 4 -V test_query.tsv
```

## Caveats

The program was made for use with SNP annotation pipelines.
Therefore, the query string only accepts a single position, not position ranges.

The code performs simple string comparison for an exact match.
