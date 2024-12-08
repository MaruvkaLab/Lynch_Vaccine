# -------------------------------------------------------------------------
# LOAD LIBRARIES
# -------------------------------------------------------------------------
library(googlesheets4)
library(ggh4x)
library(ggpubr)
library(reshape2)
library(EnvStats)
library(see)
# -------------------------------------------------------------------------
# DATA IMPORT
# -------------------------------------------------------------------------
#sig names
Indel_signatures.v3.2  <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13","ID14","ID15","ID16","ID17","ID18")
SBS_cosmic_v3.3_subset <- c("SBS1","SBS5","SBS6","SBS10a","SBS10b","SBS14","SBS15","SBS20","SBS21","SBS26")
SBS_cosmic_v3.3        <- c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d","SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14","SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33","SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40","SBS41","SBS42","SBS43","SBS44","SBS45","SBS46","SBS47","SBS48","SBS49","SBS50","SBS51","SBS52","SBS53","SBS54","SBS55","SBS56","SBS57","SBS58","SBS59","SBS60","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91","SBS92","SBS93","SBS94","SBS95")
# -------------------------------------------------------------------------
Merged <- as.data.frame(readxl::read_xlsx("supp.tables/Supplementary Table 1 - Metadata.xlsx",na = "NA"))
# select samples and columns
Merged_processed   <- Merged[!is.na(Merged$`total SNVs`), c("Lynch_germLine_mut","Histology","class.", "deconstructSigs_cosmic_v3.3.exonic", "deconstructSigs_cosmic_v3.3.subsets","Sig.ID.s COSMIC_v3.2" )]
# -------------------------------------------------------------------------
# PART 1: LIMITED NUMBER OF TESTED REFERENCE SIGNATURES (SUBSET)
# -------------------------------------------------------------------------
# Extract subset signatures
Sigs_subset <- strsplit(Merged_processed$deconstructSigs_cosmic_v3.3.subsets, split = ";")
Sigs_ID_sub <- SBS_cosmic_v3.3_subset
len_sub <- length(SBS_cosmic_v3.3_subset)
# Add signature columns to Merged_processed for the subset
for (i in seq_len(len_sub)) { Merged_processed[, Sigs_ID_sub[i]] <- as.numeric(sapply(Sigs_subset, `[`, i), 3) }
# Define used SBS signatures for subset
used_SBS_sub <- c(
  "SBS1","SBS2","SBS5","SBS6","SBS10a","SBS10b","SBS10c","SBS10d","SBS13","SBS14",
  "SBS15","SBS20","SBS21","SBS26",
  # Using indices from Sigs_ID_sub:
  Sigs_ID_sub[c(1,2,5,10,13,14,15,20,21,26)]
)
SIGS_sub <- colnames(Merged_processed)[(ncol(Merged_processed)-(len_sub-1)):ncol(Merged_processed)]
# Melt the data for subset
Merged_m_sub <- melt(Merged_processed[,c("Lynch_germLine_mut","Histology","class.",SIGS_sub)],variable.name = "sigs",  value.name = "weight")
# Prepare data for plotting for the subset scenario
selected_sub <- c("SBS6","SBS26") #select signature to plot - - 
p <- Merged_m_sub
p$Lynch_germLine_mut[is.na(p$Lynch_germLine_mut)] <- ""
p$class.[p$class. == "Lynch_like"] <- "Lynch-like"
p <- p[p$sigs %in% selected_sub, ]
# Clean up histology labels
p$Histology <- sub("CARCINOMA", "Carcinoma", p$Histology)
p$Histology <- sub("ADENOMA", "Adenoma", p$Histology)
# Clean up class labels
p$class. <- sub("Sporadic", "Sporadic-MSI", p$class.)
# Create a grouping variable
p$group <- with(p, paste0(class., "\n", Histology, "\n", Lynch_germLine_mut))
p <- p[p$group != "Lynch\nCarcinoma\n", ] #remove Lynch with Un-Identified gene
p$group <- sub("UNK","", p$group)
# Ordering the groups
p$group <- factor(p$group, levels = unique(p$group)[c(9,7,8,4,3,2,5,6,1)], ordered = TRUE)
# -------------------------------------------------------------------------
# PART 2: ALL SIGNATURES (EXONIC)
# -------------------------------------------------------------------------
selected_all <- "SBS54" #select signature to plot - - 
Sigs_all <- strsplit(Merged_processed$deconstructSigs_cosmic_v3.3.exonic, split = ";")
Sigs_ID_all <- SBS_cosmic_v3.3
len_all <- length(SBS_cosmic_v3.3)
# Add signature columns to Merged_processed for the full set
for (i in seq_len(len_all)) { Merged_processed[, Sigs_ID_all[i]] <- as.numeric(sapply(Sigs_all, `[`, i), 3)}

SIGS_all <- colnames(Merged_processed)[(ncol(Merged_processed)-(len_all-1)):ncol(Merged_processed)]

# Melt the data for full set
p2 <- melt(  Merged_processed[,c("Lynch_germLine_mut","Histology","class.",SIGS_all)],variable.name = "sigs",value.name = "weight")
# Prepare data for plotting
p2$Lynch_germLine_mut[is.na(p2$Lynch_germLine_mut)] <- ""
p2$class.[p2$class. == "Lynch_like"] <- "Lynch-like"
p2 <- p2[p2$sigs %in% selected_all, ]
# Clean up histology labels
p2$Histology <- sub("CARCINOMA", "Carcinoma", p2$Histology)
p2$Histology <- sub("ADENOMA", "Adenoma", p2$Histology)
# Clean up class labels
p2$class. <- sub("Sporadic", "Sporadic-MSI", p2$class.)
# Grouping variable
p2$group <- with(p2, paste0(class., "\n", Histology, "\n", Lynch_germLine_mut))
p2 <- p2[p2$group != "Lynch\nCarcinoma\n", ]
p2$group <- sub("UNK","", p2$group)
p2$group <- factor(p2$group, levels = unique(p2$group)[c(9,7,8,4,3,2,5,6,1)], ordered = TRUE)

# Merge the two plot data frames
p_merged <- rbind(p, p2)
p_merged$group <- factor(p_merged$group, levels = unique(p_merged$group)[c(9,7,8,4,3,2,5,6,1)], ordered = TRUE)
p_merged <- p_merged[!is.na(p_merged$group),]
# -------------------------------------------------------------------------
# PLOT - SBS signatures
# -------------------------------------------------------------------------
PLOT_mutual <- ggplot(p_merged, aes(group, weight)) +
  geom_violinhalf(
    flip = TRUE,
    alpha = 0.5,
    aes(fill = paste0(class., "\n", Histology, "\n", Lynch_germLine_mut == "MSH6"))
  ) +
  geom_jitter(
    shape = 21, size = 0.4, width = 0.15,
    aes(fill = paste0(class., "\n", Histology, "\n", Lynch_germLine_mut == "MSH6")),
    height = 0, alpha = 0.5
  ) +
  geom_boxplot(outlier.shape = NA, width = 0.18, alpha = 0.65, size = 0.3) +
  facet_wrap(~sigs, scales = "free_y", ncol = 1, strip.position = "right") +
  stat_n_text() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab("Fraction of mutations attributed to signature") +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(minor_breaks = NULL, limits = c(0, 0.5)),
      scale_y_continuous(minor_breaks = NULL, limits = c(0, 0.5)),
      scale_y_continuous(minor_breaks = NULL, limits = c(0, 0.225))
    )
  )

# -------------------------------------------------------------------------
# PLOT - Indel signatures
# -------------------------------------------------------------------------
load("ID_sig.MutationalPatterns.Rdata")
ggplot(TMP,aes(x=group,y=Weight))+
  geom_violinhalf(alpha=.25,width=1.3,flip = T,aes(fill=paste0(Class,Histology,`MMR gene`=="MSH6")))+
  geom_jitter(shape=21,size=.4,width=.1,height = 0,alpha=.5,aes(fill=paste0(Class,Histology,`MMR gene`=="MSH6")))+
  geom_boxplot(outlier.shape = NA,width=.125,alpha=.75,col="black",size=.4)+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1))+
  xlab(NULL)+
  facet_wrap(~Signature,ncol = 1,scales = "free_y",strip.position = "right")+
  ylab("Relative signature contribution")+
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)), 
               vjust=.15,hjust = -.575,color = "black",size=3)+
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, .65)),
      scale_y_continuous(limits = c(0, 0.8)),
      scale_y_continuous(limits = c(0, 0.6)),
      scale_y_continuous(limits = c(0, .5))))+
  #reduce y axis label size
  theme(axis.text.y = element_text(size = 5))+
  theme(panel.spacing = unit(0.15, "lines"))
