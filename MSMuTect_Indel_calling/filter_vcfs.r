library(data.table)

# Read VEP output and filter specific mutation types
VEP <- fread(VEP_FILE, data.table = FALSE, skip = 38)
VEP <- VEP[VEP$INFO %in% c("frameshift_variant", "inframe_deletion", "protein_altering_variant"), ]
write.table(VEP, file = paste0(OUTPUT_DIR, "/", VEP_NAME, ".filtered.vcf"), quote = FALSE, row.names = FALSE)
