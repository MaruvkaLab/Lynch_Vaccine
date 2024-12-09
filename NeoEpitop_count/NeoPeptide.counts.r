############################################################
# Author: Hagya Ladany
# Date: 9-Dec-2024
# Note: This script was modified for publication purposes.
############################################################
# Load NeoEpitope count data (calculated by MSMuTect -> Ensemble VEP -> netMHCpan4.1 via pCAVseq wrapper)
load("nature_cancer/github/Lynch_Vaccine/NeoEpitop_count/NeoEpitop_count.Rdata")
# Load metadata
MetaData <- as.data.frame(readxl::read_xlsx("nature_cancer/supp.tables/Supplementary Table 1 - Metadata.xlsx"))
# Set coverage thresholds
Cov_Thresh <- c(TUMOR = 25, NORMAL = 18)
# Select samples that meet coverage thresholds
SELECT <- which(binder_unique_per_mut$Sample %in% 
                MetaData[MetaData$mean_cov.T > Cov_Thresh["TUMOR"] & 
                         MetaData$mean_cov.N > Cov_Thresh["NORMAL"],"Paper ID"])
# Melt the data for plotting
binder_unique_per_mut.m <- reshape2::melt(binder_unique_per_mut[SELECT, ])
# Add metadata columns
binder_unique_per_mut.m <- cbind(binder_unique_per_mut.m,MetaData[ match(binder_unique_per_mut.m$Sample, MetaData$`Paper ID`),
                                                                  c("class.", "Histology", "Lynch_germLine_mut")])
# Set factor levels for germline mutation
binder_unique_per_mut.m$Lynch_germLine_mut <- factor(binder_unique_per_mut.m$Lynch_germLine_mut,
                                                     levels = c("MSH6", "MLH1", "PMS2", "MSH2", "UNK"))
# Create paper groups
binder_unique_per_mut.m$group <- with(binder_unique_per_mut.m,
                                      paste0(class., "\n", Histology, "\n", Lynch_germLine_mut))
# Standardize naming conventions in the groups
binder_unique_per_mut.m$group <- sub("CARCINOMA", "Carcinoma", binder_unique_per_mut.m$group)
binder_unique_per_mut.m$group <- sub("ADENOMA", "Adenoma", binder_unique_per_mut.m$group)
binder_unique_per_mut.m$group <- sub("Lynch_like", "Lynch-like", binder_unique_per_mut.m$group)
binder_unique_per_mut.m$group <- sub("NA", "", binder_unique_per_mut.m$group)
binder_unique_per_mut.m$group <- sub("Sporadic", "Sporadic-MSI", binder_unique_per_mut.m$group)
# Correct typo in referencing 'group' column
binder_unique_per_mut.m$group <- factor(binder_unique_per_mut.m$group,
                                        levels = c("Lynch\nAdenoma\nMSH6", 
                                                   "Lynch\nAdenoma\nMLH1", 
                                                   "Lynch\nAdenoma\nMSH2", 
                                                   "Lynch\nCarcinoma\nMSH6", 
                                                   "Lynch\nCarcinoma\nMLH1", 
                                                   "Lynch\nCarcinoma\nMSH2", 
                                                   "Lynch\nCarcinoma\nPMS2", 
                                                   "Lynch-like\nCarcinoma\n", 
                                                   "Sporadic-MSI\nCarcinoma\n"))
# Plot using ggplot2
ggplot(  subset(binder_unique_per_mut.m,!Lynch_germLine_mut %in% "UNK" &
                                        !is.na(class.) &
                                        !(class. == "Lynch" &
                                        is.na(Lynch_germLine_mut) &
                                        Histology == "CARCINOMA")),aes(x = group, y = value))+
  geom_violinhalf(alpha = .25,aes(fill = paste0(class., Histology, Lynch_germLine_mut == "MSH6")),width = 1.3,flip = TRUE)+
  geom_jitter(color = "purple4",aes(fill = paste0(class., Histology, Lynch_germLine_mut == "MSH6")),shape = 21,size = .8,width = .1,height = 0)+
  geom_boxplot(outlier.shape = NA,width = .125,alpha = .75)+
  theme_bw()+
  ylab("No. Neoepitopes")+
  labs(caption = "70nM")+
  scale_y_continuous(trans = "pseudo_log",breaks = c(0, 2, 4, 10, 33, 100, 333, 1000))+
  stat_n_text()+
  theme(legend.position = "none",panel.spacing = unit(0, 'lines'),axis.ticks.x = element_blank())+
  stat_summary(fun = median,geom = "text",aes(label = round(..y.., 2)),vjust = .15,hjust = -.3,color = "black",size = 3)+
  xlab(NULL)
