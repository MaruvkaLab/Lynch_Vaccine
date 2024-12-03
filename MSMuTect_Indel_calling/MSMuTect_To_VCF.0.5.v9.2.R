#!/usr/bin/env Rscript --vanilla

# MSMuTect_To_VCF.0.5.v9.2.R
# Converts MSMuTect output to VCF format with various filtering options.

# Log the start of the script
cat(paste0(substr(Sys.time(), 1, 19), "  @@@ MSMuTect0.5 to VCF - start:\n"))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Built-in filters:
# - Remove variants on sex chromosomes (X & Y)
# - Remove duplicated variants
# - Remove variants where the tumor allele matches the reference genome

# Display raw input arguments
cat(paste0("\tRaw input:", args, "\n\n"))
cat(paste0(substr(Sys.time(), 1, 19), "  @@@ get arguments:\n"))

# Ensure exactly 7 arguments are provided
if (length(args) != 7) {
  stop(
    paste0(
      "\n7 arguments must be supplied (you entered ", length(args), "):\n",
      "1. Input file (char)\n",
      "2. Input Path (char)\n",
      "3. Output Path (char)\n",
      "4. Single Sample Mode (Logical T/F)\n",
      "5. Filtering mode (int 1..8):\n",
      "\t1. Filter all non REF normals,\n",
      "\t2. Filter normals homozygous to ALT,\n",
      "\t3. Filter normals heterozygous to ALT,\n",
      "\t4. No filters,\n",
      "\t5. Exclude LOH: Remove variants where the normal is heterozygous and the tumor is homozygous to one of the normal alleles.\n",
      "\t6. Remove all normals with heterozygous samples.\n",
      "\t7. Remove all normals with heterozygous or homozygous ALT alleles & all tumors that are not diploid heterozygous to REF.\n",
      "\t8. Keep only tumor alleles that contain novel variants (not seen in REF or Normal).\n",
      "6. Output filtered MSMuTect file with VCF (Logical T/F)\n",
      "7. Split multi-allelic loci in VCF to single lines (Logical T/F)\n\n",
      call. = FALSE
    )
  )
}

# Assign arguments to variables
INPUT <- args[1]
inputPath <- args[2]
outputPath <- args[3]
OneSample <- as.logical(args[4])
mode <- as.integer(args[5])
output_MSMuTect <- as.logical(args[6])
singleLine <- as.logical(args[7])

# Display processed inputs
cat("Processed input:\n")
cat(paste0("\t1) Input file:.......................\t", INPUT, "\n"))
cat(paste0("\t2) Input Path:.......................\t", inputPath, "\n"))
cat(paste0("\t3) Output Path:......................\t", outputPath, "\n"))
cat(paste0("\t4) Single Sample Mode:...............\t", OneSample, "\n"))
cat(paste0("\t5) Normal filter mode:................\t", mode, "\n"))
cat(paste0("\t6) Output MSMuTect file with filters: \t", output_MSMuTect, "\n"))
cat(paste0("\t7) Separate multi-allelic variants:....\t", singleLine, "\n\n"))

# Additional arguments
remove_Sex_Chr <- TRUE

# Define filtering modes
MODES <- c(
  "1) Filter all non REF normals",
  "2) Filter normals homozygous to ALT",
  "3) Filter normals heterozygous to ALT",
  "4) No filters",
  "5) Exclude LOH: Remove variants where the normal is heterozygous and the tumor is homozygous to one of the normal alleles.",
  "6) Remove all normals with heterozygous samples.",
  "7) Remove all normals with heterozygous or homozygous ALT alleles & all tumors that are not diploid heterozygous to REF.",
  "8) Keep only tumor alleles that contain novel variants (not seen in REF or Normal)."
)

# Warning for incompatible options
if (output_MSMuTect & !singleLine) {
  warning("Can't export MSMuTect file when split multi-allelic option is on (argument 7==TRUE).")
}

### --- Import MSMuTect Input --- ###
setwd(inputPath)
options(stringsAsFactors = FALSE)
MSMuTect_output <- read.delim(INPUT)
invisible(gc())

### --- Generate VCF File --- ###
cat(paste0(substr(Sys.time(), 1, 19), "  @@@ generate VCF header\n"))

# Initialize VCF output dataframe
VCF_output <- data.frame(
  "#CHROM" = MSMuTect_output$CHROMOSOME,
  POS = MSMuTect_output$START,
  ID = ".",
  REF = NA,
  ALT = NA,
  QUAL = ".",
  FILTER = "PASS",
  INFO = paste0("LOCUS=", paste(MSMuTect_output$CHROMOSOME, MSMuTect_output$START, MSMuTect_output$END, MSMuTect_output$PATTERN, floor(MSMuTect_output$REFERENCE_REPEATS), sep = ":")),
  FORMAT = "GT",
  stringsAsFactors = FALSE
)

# Reference repeats
REFERENCE_REPEATS <- floor(MSMuTect_output$REFERENCE_REPEATS)

# Generate REF and ALT columns
cat(paste0(substr(Sys.time(), 1, 19), "  @@@ generate VCF ALT REF columns\n"))

if (OneSample) {
  # Single sample mode
  GT <- as.matrix(cbind(REFERENCE_REPEATS, MSMuTect_output[, c("ALLELE_1", "ALLELES_2", "ALLELES_3", "ALLELES_4")]))
  
  # Remove REF identical alleles
  for (i in 2:5) {
    GT[GT[, i] == GT[, 1], i] <- NA
  }
  
  # Remove rows where all alleles are NA
  RemoveMe <- !apply(GT[, 2:5], 1, function(x) all(is.na(x)))
  GT <- GT[RemoveMe, ]
  VCF_output <- VCF_output[RemoveMe, ]
  MSMuTect_output <- MSMuTect_output[RemoveMe, ]
  REFERENCE_REPEATS <- REFERENCE_REPEATS[RemoveMe]
  
  # Adjust repeats
  GT <- as.data.frame(GT)
  GT <- t(apply(GT, 1, function(x) x - (min(x, na.rm = TRUE) - 1)))
  GT_unique <- if (nrow(GT) > 1) {
    apply(GT[, 2:5], 1, function(x) unique(na.omit(x)))
  } else {
    unique(na.omit(GT[, 2:5]))
  }
  
  GT <- cbind(
    GT[, "REFERENCE_REPEATS"],
    substr(MSMuTect_output$REFERENCE_SEQUENCE, 1, nchar(MSMuTect_output$PATTERN))
  )
  colnames(GT) <- c("REF", "Base")
  
  ALT_ALLELS <- sapply(1:length(GT_unique), function(ROW) strrep(GT[ROW, "Base"], GT_unique[[ROW]]))
  VCF_output$ALT <- sapply(ALT_ALLELS, function(x) paste0(na.omit(x), collapse = ","))
  
  # Generate REF column
  VCF_output$REF <- sapply(1:nrow(GT), function(ROW) {
    substr(MSMuTect_output$REFERENCE_SEQUENCE[ROW], 1, nchar(GT[ROW, "Base"]) * as.numeric(GT[ROW, "REF"]))
  })
  
} else {
  # Regular run (Normal / Tumor)
  GT <- as.matrix(cbind(
    REFERENCE_REPEATS,
    MSMuTect_output[, c("NORMAL_ALLELE_1", "NORMAL_ALLELES_2", "NORMAL_ALLELES_3", "NORMAL_ALLELES_4",
                        "TUMOR_ALLELE_1", "TUMOR_ALLELES_2", "TUMOR_ALLELES_3", "TUMOR_ALLELES_4")]
  ))
  
  # Remove REF identical alleles
  for (i in 2:9) {
    GT[GT[, i] == GT[, 1], i] <- NA
  }
  
  GT <- as.data.frame(GT)
  GT <- t(apply(GT, 1, function(x) x - (min(x, na.rm = TRUE) - 1)))
  
  GT_unique <- if (nrow(GT) > 1) {
    apply(GT[, 2:9], 1, function(x) unique(na.omit(x)))
  } else {
    unique(na.omit(GT[, 2:9]))
  }
  
  GT <- cbind(
    GT[, "REFERENCE_REPEATS"],
    substr(MSMuTect_output$REFERENCE_SEQUENCE, 1, nchar(MSMuTect_output$PATTERN))
  )
  colnames(GT) <- c("REF", "Base")
  
  ALT_ALLELS <- sapply(1:length(GT_unique), function(ROW) strrep(GT[ROW, "Base"], GT_unique[[ROW]]))
  VCF_output$ALT <- sapply(ALT_ALLELS, function(x) paste0(na.omit(x), collapse = ","))
  
  # Generate REF column
  VCF_output$REF <- sapply(1:nrow(GT), function(ROW) {
    substr(MSMuTect_output$REFERENCE_SEQUENCE[ROW], 1, nchar(GT[ROW, "Base"]) * as.numeric(GT[ROW, "REF"]))
  })
}

# Assign genotypes
if (OneSample) {
  GT <- as.matrix(cbind(REFERENCE_REPEATS, MSMuTect_output[, c("ALLELE_1", "ALLELES_2", "ALLELES_3", "ALLELES_4")]))
  GT_unique <- apply(GT[, 1:5], 1, function(x) unique(na.omit(x)))
  GT <- GT[, -1]
  GT <- as.data.frame(GT)
  GT_score <- t(sapply(1:nrow(GT), function(x) match(GT[x, ], GT_unique[[x]]) - 1))
  colnames(GT_score) <- colnames(GT)
  
  # Fill missing alleles assuming diploid
  GT_score[is.na(GT_score[, 2]), 2] <- GT_score[is.na(GT_score[, 2]), 1]
  
  GT_Normal <- paste(GT_score[, 1], GT_score[, 2], GT_score[, 3], GT_score[, 4], sep = "/")
  GT_Normal <- gsub("/NA", "", GT_Normal)
  VCF_output$TUMOR <- GT_Normal
  
} else {
  GT <- as.matrix(cbind(
    REFERENCE_REPEATS,
    MSMuTect_output[, c("NORMAL_ALLELE_1", "NORMAL_ALLELES_2", "NORMAL_ALLELES_3", "NORMAL_ALLELES_4",
                        "TUMOR_ALLELE_1", "TUMOR_ALLELES_2", "TUMOR_ALLELES_3", "TUMOR_ALLELES_4")]
  ))
  
  GT_unique <- apply(GT[, 1:9], 1, function(x) unique(na.omit(x)))
  GT <- GT[, -1]
  GT <- as.data.frame(GT)
  GT_score <- t(sapply(1:nrow(GT), function(x) match(GT[x, ], GT_unique[[x]]) - 1))
  colnames(GT_score) <- colnames(GT)
  
  # Fill missing alleles assuming diploid
  GT_score[is.na(GT_score[, 2]), 2] <- GT_score[is.na(GT_score[, 2]), 1]
  GT_Normal <- paste(GT_score[, 1], GT_score[, 2], GT_score[, 3], GT_score[, 4], sep = "/")
  GT_Normal <- gsub("/NA", "", GT_Normal)
  
  GT_score[is.na(GT_score[, 6]), 6] <- GT_score[is.na(GT_score[, 6]), 5]
  GT_Tumor <- paste(GT_score[, 5], GT_score[, 6], GT_score[, 7], GT_score[, 8], sep = "/")
  GT_Tumor <- gsub("/NA", "", GT_Tumor)
  GT_Tumor[nchar(GT_Tumor) == 1] <- paste0(GT_Tumor[nchar(GT_Tumor) == 1], "/", GT_Tumor[nchar(GT_Tumor) == 1])
  
  VCF_output$NORMAL <- GT_Normal
  VCF_output$TUMOR <- GT_Tumor
}

# Remove sex chromosomes if specified
if (remove_Sex_Chr) {
  if (output_MSMuTect & !singleLine) {
    MSMuTect_output <- MSMuTect_output[!(VCF_output$`#CHROM` %in% c("X", "Y")), ]
  }
  VCF_output <- VCF_output[!(VCF_output$`#CHROM` %in% c("X", "Y")), ]
}

# Order VCF output
if (output_MSMuTect & !singleLine) {
  MSMuTect_output <- MSMuTect_output[order(as.numeric(VCF_output$`#CHROM`), as.numeric(VCF_output$POS)), ]
  MSMuTect_output <- MSMuTect_output[!duplicated(VCF_output[, 1:5]), ]  # Export unique variants only
}
VCF_output <- VCF_output[order(as.numeric(VCF_output$`#CHROM`), as.numeric(VCF_output$POS)), ]
VCF_output_unique <- VCF_output[!duplicated(VCF_output[, 1:5]), ]  # Export unique variants only

setwd(outputPath)

# Function to remove tumor alleles identical to REF
RM_homoREF <- function() {
  MSMuTect_output <<- MSMuTect_output[VCF_output_unique$TUMOR != "0/0", ]
  VCF_output_unique <<- VCF_output_unique[VCF_output_unique$TUMOR != "0/0", ]
}

# Function to split multi-allelic loci into single alleles using bcftools
singleLines <- function() {
  if (singleLine) {
    NAME <- paste0(INPUT, '_', Sys.Date(), '_F', mode, '.vcf')
    cat(paste0(substr(Sys.time(), 1, 19), "  @@@ split multi-allelic loci to single lines... bgzip..."))
    
    system(paste0("bgzip ", outputPath, "/", DIRn, "/", NAME))
    cat(" tabix...")
    system(paste0("tabix ", outputPath, "/", DIRn, "/", NAME, ".gz"))
    cat(" split calls...\n")
    
    system(paste0("bcftools norm -m -any ", outputPath, "/", DIRn, "/", NAME, ".gz > ",
                  outputPath, "/", DIRn, "/", "SingleLine_", NAME))
    
    cat("\t\t\t\t\t\t\t clean intermediate files...\n")
    file.remove(paste0(outputPath, "/", DIRn, "/", NAME, ".gz"))
    
    cat("\t\t\t\t\t\t\t read file again to modify tumor and normal columns\n")
    FNAME <- paste0(outputPath, "/", DIRn, "/", "SingleLine_", NAME)
    bcfOUT <- read.delim(FNAME, comment = "#", header = FALSE)
    colnames(bcfOUT) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
    
    bcfOUT$TUMOR <- "0/1"
    if (sum(grepl("1/0", bcfOUT$NORMAL)) > 1) {
      bcfOUT[grepl("1/0", bcfOUT$NORMAL), "NORMAL"] <- "0/1"
    }
    
    system(paste0("grep '##' ", FNAME, " > ", sub("SingleLine_", "Sl_", FNAME)))
    
    suppressWarnings(
      write.table(bcfOUT, sub("SingleLine_", "Sl_", FNAME), sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = TRUE, append = TRUE)
    )
    
    file.remove(FNAME)                                           # Clean intermediate files
    file.remove(paste0(outputPath, "/", DIRn, "/", NAME, ".gz.tbi"))  # Clean index file
    
    cat("\t\t\t\t\t\t split multi-allelic loci - Finished\n")
  }
}

# Function to save VCF and MSMuTect files
Save_Files <- function() {
  cat(paste0("\n", substr(Sys.time(), 1, 19), "  @@@ Saving VCF file..."))
  
  # Write VCF header
  write.table(
    VCF_header,
    paste0(outputPath, "/", DIRn, "/", INPUT, '_', Sys.Date(), '_F', mode, '.vcf'),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  # Ensure all columns are character type
  for (col in 1:ncol(VCF_output_unique)) {
    class(VCF_output_unique[, col]) <- "character"
  }
  
  # Write VCF data
  write.table(
    VCF_output_unique,
    paste0(outputPath, "/", DIRn, "/", INPUT, '_', Sys.Date(), '_F', mode, '.vcf'),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE
  )
  
  # Validate VCF file creation
  if (file.exists(paste0(outputPath, "/", DIRn, "/", INPUT, '_', Sys.Date(), '_F', mode, '.vcf'))) {
    cat(paste0(substr(Sys.time(), 1, 19), "  VCF file created.\n"))
  } else {
    stop("Could not validate VCF file on disk. Check write permissions.")
  }
  
  # Save filtered MSMuTect file if required
  if (output_MSMuTect & !singleLine) {
    cat(paste0(substr(Sys.time(), 1, 19), "  @@@ Saving filtered MSMuTect file...\n"))
    write.table(
      MSMuTect_output,
      file = paste0(outputPath, "/", DIRn, "/", gsub(pattern = "full.mut", replacement = "filtered", x = INPUT)),
      quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
    )
    
    # Validate MSMuTect file creation
    if (file.exists(paste0(outputPath, "/", DIRn, "/", gsub(pattern = "full.mut", replacement = "filtered", x = INPUT)))) {
      cat(paste0(substr(Sys.time(), 1, 19), "  MSMuTect.filtered file created.\n"))
    } else {
      stop("Something is wrong, file not found?!")
    }
  }
}

### --- Filtering Based on Mode --- ###
if (!OneSample) {
  cat(paste0(substr(Sys.time(), 1, 19), "  @@@ Normal filtering step: "))
  
  if (mode == 1) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine) {
      MSMuTect_output <- MSMuTect_output[VCF_output_unique$NORMAL %in% c("0/0", "0|0"), ]
    }
    VCF_output_unique <- VCF_output_unique[VCF_output_unique$NORMAL %in% c("0/0", "0|0"), ]
    DIRn <- "Keep_Only_Homozygous_to_REF"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 2) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine) {
      MSMuTect_output <- MSMuTect_output[grepl("0", VCF_output_unique$NORMAL), ]
    }
    VCF_output_unique <- VCF_output_unique[grepl("0", VCF_output_unique$NORMAL), ]
    DIRn <- "Filter_normals_homozygous_to_ALT"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 3) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine) {
      MSMuTect_output <- MSMuTect_output[!grepl(pattern = "0/[1-9]|[1-9]/0", VCF_output_unique$NORMAL), ]
    }
    VCF_output_unique <- VCF_output_unique[!grepl(pattern = "0/[1-9]|[1-9]/0", VCF_output_unique$NORMAL), ]
    DIRn <- "Filter_normals_heterozygous_to_ALT"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 4) {
    DIRn <- "N_Not_Filtered"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 5) {
    RM_homoREF()
    RemoveMe <- -as.numeric(names(which(
      colSums(apply(
        VCF_output_unique[!VCF_output_unique$NORMAL %in% c("0/0", "1/1") &
                           nchar(VCF_output_unique$TUMOR) < 4,
                         c("NORMAL", "TUMOR")],
        1,
        function(Row) any(strsplit(Row[2], "/")[[1]] %in% strsplit(Row[1], "/")[[1]])
      )) == 2
    )))
    VCF_output_unique <- VCF_output_unique[RemoveMe, ]
    if (output_MSMuTect & !singleLine) {
      MSMuTect_output <- MSMuTect_output[RemoveMe, ]
    }
    DIRn <- "Exclude_LOH"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 6) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine) {
      MSMuTect_output <- MSMuTect_output[sapply(strsplit(VCF_output_unique$NORMAL, "/"), function(X) X[1] == X[2]), ]
    }
    VCF_output_unique <- VCF_output_unique[sapply(strsplit(VCF_output_unique$NORMAL, "/"), function(X) X[1] == X[2]), ]
    DIRn <- "Exclude_normal_heterozygous"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 7) {
    RM_homoREF()
    if (output_MSMuTect & !singleLine) {
      MSMuTect_output <- MSMuTect_output[VCF_output_unique$NORMAL == "0/0" &
                                         VCF_output_unique$TUMOR %in% c("0/1", "1/0"), ]
    }
    VCF_output_unique <- VCF_output_unique[VCF_output_unique$NORMAL == "0/0" &
                                           VCF_output_unique$TUMOR %in% c("0/1", "1/0"), ]
    DIRn <- "Keep_N_Homozygous_to_REF_&_T_diploid_heterozygous_tumors"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else if (mode == 8) {
    RM_homoREF()
    RemoveMe <- -as.numeric(names(which(
      colSums(apply(
        VCF_output_unique[!VCF_output_unique$NORMAL %in% c("0/0", "1/1") &
                           nchar(VCF_output_unique$TUMOR) < 4,
                         c("NORMAL", "TUMOR")],
        1,
        function(Row) any(strsplit(Row[2], "/")[[1]] %in% strsplit(Row[1], "/")[[1]])
      )) == 2
    )))
    
    if (length(RemoveMe) > 1) {
      VCF_output_unique <- VCF_output_unique[RemoveMe, ]
      MSMuTect_output <- MSMuTect_output[RemoveMe, ]
    }
    
    RemoveMe <- sapply(1:nrow(VCF_output_unique), function(ROW) {
      all(unlist(strsplit(VCF_output_unique$TUMOR[ROW], "/")) %in% 
            c(unlist(strsplit(VCF_output_unique$NORMAL[ROW], "/")), "0"))
    })
    
    VCF_output_unique <- VCF_output_unique[!RemoveMe, ]
    MSMuTect_output <- MSMuTect_output[!RemoveMe, ]
    
    DIRn <- "novel_only"
    dir.create(file.path(DIRn), showWarnings = FALSE)
    cat(MODES[mode])
    Save_Files()
    singleLines()
    
  } else {
    stop("\nBad input in 'mode' argument. Valid values: 1..8.\n")
  }
}

# Completion message
cat("\n\n\tDone. Shalom aleichem.\n\n")

### --- Generate VCF Header --- ###
if (OneSample) {
  VCF_header <- paste0(
    '##fileformat=VCFv4.2
##fileDate=', Sys.Date(), '
##INFO=<ID=LOCUS,Number=.,Type=String,Description="Microsatelite locus taken from MSMuTect output">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##source=Yosef Maruvka`s lab:www.MaruvkaLab.com
##Tool=github.com/Hagay-Ladany/MSMuTect-0.5-to-VCF Input:', 
    INPUT, ',', inputPath, ',', outputPath, ',', OneSample, ',', mode, ',', output_MSMuTect, ',', singleLine, '
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR'
  )
} else {
  VCF_header <- paste0(
    '##fileformat=VCFv4.2
##fileDate=', Sys.Date(), '
##INFO=<ID=LOCUS,Number=.,Type=String,Description="Microsatelite locus taken from MSMuTect output">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##Tool=github.com/Hagay-Ladany/MSMuTect-0.5-to-VCF Input:', 
    INPUT, ',', inputPath, ',', outputPath, ',', OneSample, ',', mode, ',', output_MSMuTect, ',', singleLine, '
##source=Yosef Maruvka`s lab:www.MaruvkaLab.com
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR'
  )
}
