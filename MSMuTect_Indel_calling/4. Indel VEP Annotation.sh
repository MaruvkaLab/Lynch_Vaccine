#!/bin/bash

export PATH=/path/to/anaconda3/bin:$PATH
# Run VEP annotation
/path/to/vep \
 --cache \
 --assembly GRCh37 \
 --dir_cache /path/to/vep_cache/ \
 --dir_plugins /path/to/plugins/ \
 --input_file $FILE_NAME.gz \
 --output_file /output/path/VEP_output/$FILE_NAME.VEP.vcf \
 --plugin Frameshift \
 --plugin Wildtype \
 --plugin Downstream \
 --plugin NMD \
 --format vcf \
 --terms SO \
 --tsl \
 --hgvs \
 --offline \
 --fasta /path/to/reference.fasta \
 --vcf \
 --symbol \
 --pick
