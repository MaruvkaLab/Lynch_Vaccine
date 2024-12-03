#!/bin/bash
# Input: Sample information (Normal and Tumor BAM files)
IFS=";" read NORMAL TUMOR <<< $SAMPLE
# Set environment paths
export PATH=/path/to/anaconda3/bin:$PATH
export PYTHONPATH="/path/to/MSMuTect"
echo "Normal: $NORMAL.bam"
echo "Tumor: $TUMOR.bam"
# Run MSMuTect
python3 /path/to/MSMuTect/src/Entry/main.py \
 -l /path/to/phobos_results/all.sorted.fa.phobos \
 -T $TUMOR.bam \
 -N $NORMAL.bam \
 -O /output/path/MSMuTect_output/$NORMAL.AND.$TUMOR \
 -m -A -H -c 10
echo "MSMuTect finished"
