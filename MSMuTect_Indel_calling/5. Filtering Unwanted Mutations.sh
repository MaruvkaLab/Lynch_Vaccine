#!/bin/bash

export PATH=/path/to/R/bin:$PATH

# Run filtering script
Rscript --vanilla filter_vcfs.r $VEP_FILE $VEP_NAME $OUTPUT
