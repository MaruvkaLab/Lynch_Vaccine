#!/bin/bash

# Filter Indel mutations with status "M"
cat $SAMPLE.tsv | head -n 1 > Called.only/$SAMPLE.MutOnly
awk -F'\t' '$49 == "M"' $SAMPLE.tsv >> Called.only/$SAMPLE.MutOnly
