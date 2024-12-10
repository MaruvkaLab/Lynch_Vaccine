# Name: Hagya Ladany
# Date: 9-Dec-2024
# This script visualizes the evaluation of vaccine success using data from netMHCpan4.1 binding affinities
# and comprehensive scores for both off-the-shelf and HLA-based (personalized) vaccine approaches.
# this script was modified for publication purposes

# load libraries
library(ggplot2)
library(cowplot)
library(ggpubr)
library(data.table)
setwd("C:\\Users\\Bad-G\\OneDrive - Technion\\Third_Degree\\lynch vaccine\\MMRd_paper_maruvka\\nature_cancer\\github\\Lynch_Vaccine\\vaccine_visual_evaluation")

# Visualization based on netMHCpan4.1 binding affinities (off-the-shelf and HLA-based approaches)
load("Rez.noDupSAM.byBinding.score.13Nov2024.Rdata")
# select of the shelf data
General <- do.call(rbind, rez[grep("^UNI", names(rez))])
General$Group[General$Group == "ALL"] <- paste0("All n=", General$N_pat[General$Group == "ALL"][1])
General$Group[General$Group == "sporadic"] <- paste0("Sporadic n=", General$N_pat[General$Group == "sporadic"][1])
General <- General[General$Group != "Lynch.carcinoma", ]
General$Group[General$Group == "Lynch.carcinoma.nonMSH6"] <- paste0("Lynch Carcinoma\n(non-MSH6) n=", General$N_pat[General$Group == "Lynch.carcinoma.nonMSH6"][1])
General$Group[General$Group == "ADENOMA"] <- paste0("Adenomas n=", General$N_pat[General$Group == "ADENOMA"][1])
General$Group[General$Group == "MSH6.carcinoma"] <- paste0("MSH6 Carcinoma n=", General$N_pat[General$Group == "MSH6.carcinoma"][1])

SUMM <- General[General$No.Hits %in% c("1 hit", "6 hits") & General$pepGiven %in% c(10, 50, 200), ]
SUMM$Frac.pat <- round(SUMM$Frac.pat, digits = 2)

# Plot off-the-shelf approach success rates for at least one hit
F1 <- ggplot(data = General[which(General$No.Hits == "1 hit"), ], aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_point() + geom_line() + theme_classic() + theme(legend.position = "right") +
  xlab("No. peptides") + ylab("Fraction of patients") + theme(legend.title = element_blank())

# Plot off-the-shelf approach success rates for all hits
F2 <- ggplot(data = General[General$No.Hits %in% c("1 hit", "6 hits"), ], aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  theme_bw() + theme(legend.position = "right") + xlab("No. peptides") +
  facet_wrap(~No.Hits, ncol = 6) + ylab("Fraction of patients") +
  scale_y_continuous(breaks = (seq(0, 10, 2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  # Increase space between legend labels
  theme(legend.text = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(legend.key.size = unit(1, "cm")) + theme(legend.spacing.x = unit(1, 'cm')) +
  # Remove space under the plot
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
  # Remove space between plot and legend
  theme(legend.margin = margin(0, 0, 0, 0, "cm"))

my_legend <- get_legend(F2)
ggsave(as_ggplot(my_legend), filename = "General_Binding.score/General_all_hits_legend.svg", width = 15, height = 3)
ggsave(plot = F2 + theme(legend.position = "none"), filename = paste0("General_Binding.score/General_all_hits.big.", format(Sys.time(), "%d%b%Y"), ".svg"), width = 4, height = 2)

Personalized <- do.call(rbind, rez[grep("^P", names(rez))])
Personalized$Group[Personalized$Group == "ALL"] <- paste0("All n=", Personalized$N_pat[Personalized$Group == "ALL"][1])
Personalized$Group[Personalized$Group == "sporadic"] <- paste0("Sporadic n=", Personalized$N_pat[Personalized$Group == "sporadic"][1])
Personalized <- Personalized[Personalized$Group != "Lynch.carcinoma", ]
Personalized$Group[Personalized$Group == "ADENOMA"] <- paste0("Adenoma n=", Personalized$N_pat[Personalized$Group == "ADENOMA"][1])
Personalized$Group[Personalized$Group == "MSH6.carcinoma"] <- paste0("MSH6 Carcinoma n=", Personalized$N_pat[Personalized$Group == "MSH6.carcinoma"][1])
Personalized$Group[Personalized$Group == "Lynch.carcinoma.nonMSH6"] <- paste0("Lynch Carcinoma\n(non-MSH6) n=", Personalized$N_pat[Personalized$Group == "Lynch.carcinoma.nonMSH6"][1])

SUMM.personalized <- Personalized[Personalized$No.Hits %in% c("1 hit", "6 hits") & Personalized$pepGiven %in% c(10, 50, 200), ]
SUMM.personalized$Frac.pat <- round(SUMM.personalized$Frac.pat, digits = 2)

# Plot personalized approach success rates for one hit
ggplot(data = Personalized[which(Personalized$No.Hits == "1 hit"), ], aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_point() + geom_line() + theme_classic() + theme(legend.position = "right") +
  xlab("No. peptides") + ylab("Fraction of patients")

F4 <- ggplot(data = Personalized, aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_line() + theme_bw() + theme(legend.position = "right") + xlab("No. peptides") +
  facet_wrap(~No.Hits, ncol = 6) + ylab("Fraction of patients") +
  scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL)

Personalized$type <- "HLA-based"
General$type <- "Off-the-shelf"
mutualCols <- intersect(colnames(General), colnames(Personalized))
Merged.rez <- rbind(General[, mutualCols], Personalized[, mutualCols])
Merged.rez$Group <- sub("Adenomas", "Adenoma", Merged.rez$Group)
Merged.rez$Group <- sub("n=", "\nn=", Merged.rez$Group)

# Compare off-the-shelf vs. personalized for 1 hit and 6 hits
F5 <- ggplot(data = Merged.rez[Merged.rez$No.Hits %in% c("1 hit", "6 hits"), ], aes(x = pepGiven, y = Frac.pat, colour = type)) +
  geom_line(size = .85) + theme_bw() + theme(legend.position = "right") +
  xlab("No. peptides") + facet_grid(Group ~ No.Hits) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  ylab("Fraction of patients") + scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  scale_color_manual(values = c("HLA-based"="#4059AD","Off-the-shelf"="#F4B902")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(plot.margin = margin(0.1,0.1,0,0,"cm")) +
  theme(legend.margin = margin(-.25,0,0,0,"cm")) +
  theme(legend.text = element_text(margin = margin(t=0,r=10,b=0,l=0))) +
  theme(legend.key.size = unit(.8,"cm")) +
  theme(legend.spacing.x = unit(0,"cm"))

ggsave(plot = F5, filename = paste0("General_Binding.score/General.VS.Pers_all_hits.big.", format(Sys.time(), "%d%b%Y"), ".svg"), width = 4.5, height = 6)
ggsave(plot = F5, filename = paste0("General_Binding.score/General.VS.Pers_all_hits.big.", format(Sys.time(), "%d%b%Y"), ".png"), width = 4.5, height = 6)

F6 <- ggplot(data = Merged.rez, aes(x = pepGiven, y = Frac.pat, colour = type)) +
  geom_line(size = .85) + theme_bw() + theme(legend.position = "right") +
  xlab("No. peptides") + facet_grid(Group ~ No.Hits) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  ylab("Fraction of patients") + scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  scale_color_manual(values = c("HLA-based"="#4059AD","Off-the-shelf"="#F4B902")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(plot.margin = margin(0.1,0.1,0.1,0.1,"cm")) +
  theme(legend.margin = margin(-.2,0,0,0,"cm")) +
  theme(legend.text = element_text(margin = margin(t=0,r=10,b=0,l=0))) +
  theme(legend.key.size = unit(.8,"cm")) +
  theme(legend.spacing.x = unit(1,"cm"))

ggsave(plot = F6, filename = paste0("General_Binding.score/General.VS.Pers_All_hits", format(Sys.time(), "%d%b%Y"), ".svg"), width = 8.5, height = 6.2)
ggsave(plot = F6, filename = paste0("General_Binding.score/General.VS.Pers_All_hits", format(Sys.time(), "%d%b%Y"), ".png"), width = 8.5, height = 6.2)

Tmp <- Merged.rez[Merged.rez$No.Hits %in% c("1 hit", "6 hits") & Merged.rez$pepGiven %in% c(10,50,200), ]
# Reshape Merged.rez using dcast with Group as columns
Merged.rez.d <- data.table::dcast(Tmp, No.Hits + Group ~ type + pepGiven, value.var = "Frac.pat")
Merged.rez.d <- Merged.rez.d[order(Merged.rez.d$Group), ]
Merged.rez.d <- Merged.rez.d[, c(1,2,6,3,7,4,8,5)]
write.csv(Merged.rez.d, "General_Binding.score/General.VS.Pers_all_hits.csv", row.names = F)

MeanMedGeneral <- do.call(rbind, rez[grepl("^M", names(rez)) & !grepl("Personalized", names(rez))])
MeanMedGeneral$group <- sub(".Binding.scor.*", "", sub("MeanMeds.UNI.SAM.", "", rownames(MeanMedGeneral)))
rownames(MeanMedGeneral) <- NULL
MeanMedGeneral$group[MeanMedGeneral$group == "ALL"] <- grep("ALL", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedGeneral$group[MeanMedGeneral$group == "sporadic"] <- grep("sporadic", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedGeneral <- MeanMedGeneral[MeanMedGeneral$group != "Lynch.carcinoma", ]
MeanMedGeneral$group[MeanMedGeneral$group == "ADENOMA"] <- grep("ADENOMA", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedGeneral$group[MeanMedGeneral$group == "MSH6.carcinoma"] <- grep("MSH6.carcinoma", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedGeneral$group[MeanMedGeneral$group == "Lynch.carcinoma.nonMSH6"] <- grep("Lynch", ignore.case = T, Personalized$Group, value = T)[1]

ggplot(data = MeanMedGeneral, aes(x = pepGiven, y = value, colour = group)) +
  geom_line(size = .85, alpha = .7) + theme_bw() + facet_wrap(~variable, ncol=2) +
  ylab("No. hits") + xlab("No. peptides") +
  theme(plot.margin = margin(0.1,0.1,0.1,0.1,"cm"), legend.position = "none")
ggsave("General_Binding.score/General_MeanMed_.big.svg", width = 4, height = 2)

MeanMedPesonalized <- do.call(rbind, rez[grepl("^M", names(rez)) & grepl("Personalized", names(rez))])
MeanMedPesonalized$group <- sub(".BindingScor.*", "", sub("MeansMEDs.Personalized.UNI.SAM.", "", rownames(MeanMedPesonalized)))
rownames(MeanMedPesonalized) <- NULL
MeanMedPesonalized$group[MeanMedPesonalized$group == "ALL"] <- grep("ALL", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized$group[MeanMedPesonalized$group == "sporadic"] <- grep("sporadic", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized <- MeanMedPesonalized[MeanMedPesonalized$group != "Lynch.carcinoma", ]
MeanMedPesonalized$group[MeanMedPesonalized$group == "ADENOMA"] <- grep("ADENOMA", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized$group[MeanMedPesonalized$group == "MSH6.carcinoma"] <- grep("MSH6.carcinoma", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized$group[MeanMedPesonalized$group == "Lynch.carcinoma.nonMSH6"] <- grep("Lynch", ignore.case = T, Personalized$Group, value = T)[1]

ggplot(data = MeanMedPesonalized, aes(x = GivenPep, y = value, colour = group)) +
  geom_point(size = .7) + geom_line() + theme_bw() + facet_wrap(~Type) +
  ylab("No. hits") + xlab("No. peptides")
ggsave("Personalized_Binding.score/Personalized_MeanMed.svg", width = 8, height = 5)

# -------------------------------------------
# Visualization based on comprehensive score
# -------------------------------------------
load("Rez.noDupSAM.byY.score1.13Nov2024.Rdata")

General <- do.call(rbind, rez[grep("^UNI", names(rez))])
General$Group[General$Group == "ALL"] <- paste0("All n=", General$N_pat[General$Group == "ALL"][1])
General$Group[General$Group == "sporadic"] <- paste0("Sporadic n=", General$N_pat[General$Group == "sporadic"][1])
General <- General[General$Group != "Lynch.carcinoma", ]
General$Group[General$Group == "Lynch.carcinoma.nonMSH6"] <- paste0("Lynch Carcinoma\n(non-MSH6) n=", General$N_pat[General$Group == "Lynch.carcinoma.nonMSH6"][1])
General$Group[General$Group == "ADENOMA"] <- paste0("Adenomas n=", General$N_pat[General$Group == "ADENOMA"][1])
General$Group[General$Group == "MSH6.carcinoma"] <- paste0("MSH6 Carcinoma n=", General$N_pat[General$Group == "MSH6.carcinoma"][1])

SUMM <- General[General$No.Hits %in% c("1 hit", "6 hits") & General$pepGiven %in% c(10,50,200), ]
SUMM$Frac.pat <- round(SUMM$Frac.pat, digits = 2)

# Plot off-the-shelf (comprehensive score) approach success rates for one hit
F1 <- ggplot(data = General[which(General$No.Hits == "1 hit"), ], aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_point() + geom_line() + theme_classic() + theme(legend.position = "right") +
  xlab("No. peptides") + ylab("Fraction of patients") +
  theme(legend.title = element_blank())
# plot 1/6 hits
F2 <- ggplot(data = General[General$No.Hits %in% c("1 hit","6 hits"), ], aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_line(size = .85, alpha = .7) + theme_bw() + theme(legend.position = "right") +
  xlab("No. peptides") + facet_wrap(~No.Hits, ncol = 6) +
  ylab("Fraction of patients") + scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) + theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(legend.text = element_text(margin = margin(t=0,r=10,b=0,l=0))) +
  theme(legend.key.size = unit(1,"cm")) + theme(legend.spacing.x = unit(1,"cm")) +
  theme(plot.margin = margin(0.1,0.1,0.1,0.1,"cm")) + theme(legend.margin = margin(0,0,0,0,"cm"))

ggsave(plot = F2+theme(legend.position = "none"), filename = paste0("General_Y.score1/General_all_hits_legend", format(Sys.time(), "%d%b%Y"), ".svg"), width = 5, height = 3)
ggsave(as_ggplot(get_legend(F2)), filename = "General_Y.score1/General_all_hits_legend.svg", width = 15, height = 3)
# personalized approach
Personalized <- do.call(rbind, rez[grep("^P", names(rez))])
Personalized$Group[Personalized$Group == "ALL"] <- paste0("All n=", Personalized$N_pat[Personalized$Group == "ALL"][1])
Personalized$Group[Personalized$Group == "sporadic"] <- paste0("Sporadic n=", Personalized$N_pat[Personalized$Group == "sporadic"][1])
Personalized <- Personalized[Personalized$Group != "Lynch.carcinoma", ]
Personalized$Group[Personalized$Group == "ADENOMA"] <- paste0("Adenoma n=", Personalized$N_pat[Personalized$Group == "ADENOMA"][1])
Personalized$Group[Personalized$Group == "MSH6.carcinoma"] <- paste0("MSH6 Carcinoma n=", Personalized$N_pat[Personalized$Group == "MSH6.carcinoma"][1])
Personalized$Group[Personalized$Group == "Lynch.carcinoma.nonMSH6"] <- paste0("Lynch Carcinoma\n(non-MSH6) n=", Personalized$N_pat[Personalized$Group=="Lynch.carcinoma.nonMSH6"][1])

SUMM.personalized <- Personalized[Personalized$No.Hits %in% c("1 hit", "6 hits") & Personalized$pepGiven %in% c(10,50,200), ]
SUMM.personalized$Frac.pat <- round(SUMM.personalized$Frac.pat, digits = 2)

ggplot(data = Personalized[which(Personalized$No.Hits == "1 hit"), ], aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_point() + geom_line() + theme_classic() + theme(legend.position = "right") +
  xlab("No. peptides") + ylab("Fraction of patients")

F4 <- ggplot(data = Personalized, aes(x = pepGiven, y = Frac.pat, colour = Group)) +
  geom_line() + theme_bw() + theme(legend.position = "right") + xlab("No. peptides") +
  facet_wrap(~No.Hits, ncol = 6) + ylab("Fraction of patients") +
  scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL)

Personalized$type <- "HLA-based"
General$type <- "Off-the-shelf"
mutualCols <- intersect(colnames(General), colnames(Personalized))
Merged.rez <- rbind(General[, mutualCols], Personalized[, mutualCols])
Merged.rez$Group <- sub("Adenomas", "Adenoma", Merged.rez$Group)

# Compare off-the-shelf vs. personalized using comprehensive score
F5 <- ggplot(data = Merged.rez[Merged.rez$No.Hits %in% c("1 hit", "6 hits"), ], aes(x = pepGiven, y = Frac.pat, colour = type)) +
  geom_line(size = .85) + theme_bw() + theme(legend.position = "right") +
  xlab("No. peptides") + facet_grid(Group~No.Hits) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  ylab("Fraction of patients") + scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  scale_color_manual(values = c("HLA-based"="#4059AD","Off-the-shelf"="#F4B902")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(plot.margin = margin(0.1,0.1,0,0,"cm")) +
  theme(legend.margin = margin(-.25,0,0,0,"cm")) +
  theme(legend.text = element_text(margin = margin(t=0,r=10,b=0,l=0))) +
  theme(legend.key.size = unit(.8,"cm")) +
  theme(legend.spacing.x = unit(0,"cm"))

ggsave(plot = F5, filename = paste0("General_Y.score1/General.VS.Pers_1_6_hits", format(Sys.time(), "%d%b%Y"), ".svg"), width = 4.5, height = 6)
ggsave(plot = F5, filename = paste0("General_Y.score1/General.VS.Pers_1_6_hits", format(Sys.time(), "%d%b%Y"), ".png"), width = 4.5, height = 6)

Tmp <- Merged.rez[Merged.rez$No.Hits %in% c("1 hit", "6 hits") & Merged.rez$pepGiven %in% c(10,50,200), ]

F6 <- ggplot(data = Merged.rez, aes(x = pepGiven, y = Frac.pat, colour = type)) +
  geom_line(size = .85) + theme_bw() + theme(legend.position = "right") +
  xlab("No. peptides") + facet_grid(Group~No.Hits) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  ylab("Fraction of patients") + scale_y_continuous(breaks = (seq(0,10,2)/10), minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL) +
  scale_color_manual(values = c("HLA-based"="#4059AD","Off-the-shelf"="#F4B902")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(plot.margin = margin(0.1,0.1,0,0,"cm")) +
  theme(legend.margin = margin(-.25,0,0,0,"cm")) +
  theme(legend.text = element_text(margin = margin(t=0,r=10,b=0,l=0))) +
  theme(legend.key.size = unit(.8,"cm")) +
  theme(legend.spacing.x = unit(0,"cm"))

ggsave(plot = F6, filename = paste0("General_Y.score1/General.VS.Pers_All_hits", format(Sys.time(), "%d%b%Y"), ".svg"), width = 8.5, height = 6)
ggsave(plot = F6, filename = paste0("General_Y.score1/General.VS.Pers_All_hits", format(Sys.time(), "%d%b%Y"), ".png"), width = 8.5, height = 6)

Merged.rez.d <- data.table::dcast(Tmp, No.Hits + Group ~ type + pepGiven, value.var = "Frac.pat")
Merged.rez.d$Group <- factor(Merged.rez.d$Group, levels = c("All n=209","Sporadic n=101","Lynch Carcinoma\n(non-MSH6) n=49","MSH6 Carcinoma n=22","Adenoma n=14"))
Merged.rez.d <- Merged.rez.d[order(Merged.rez.d$Group), ]
Merged.rez.d <- Merged.rez.d[, c(1,2,6,3,7,4,8,5)]
write.csv(Merged.rez.d, "General_Y.score1/General.VS.Pers_all_hits_YSCORE.csv", row.names = F)

MeanMedPesonalized <- do.call(rbind, rez[grepl("^M", names(rez)) & grepl("Personalized", names(rez))])
MeanMedPesonalized$group <- sub(".Y.score..*", "", sub(".*.UNI.SAM.", "", rownames(MeanMedPesonalized)))
rownames(MeanMedPesonalized) <- NULL
MeanMedPesonalized$group[MeanMedPesonalized$group=="ALL"] <- grep("ALL", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized$group[MeanMedPesonalized$group=="sporadic"] <- grep("sporadic", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized <- MeanMedPesonalized[MeanMedPesonalized$group != "Lynch.carcinoma", ]
MeanMedPesonalized$group[MeanMedPesonalized$group=="ADENOMA"] <- grep("ADENOMA", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized$group[MeanMedPesonalized$group=="MSH6.carcinoma"] <- grep("MSH6.carcinoma", ignore.case = T, Personalized$Group, value = T)[1]
MeanMedPesonalized$group[MeanMedPesonalized$group=="Lynch.carcinoma.nonMSH6"] <- grep("Lynch", ignore.case = T, Personalized$Group, value = T)[1]

ggplot(data = MeanMedPesonalized, aes(x = GivenPep, y = value, colour = group)) +
  geom_point(size = .7) + geom_line() + theme_bw() + facet_wrap(~Type) +
  ylab("No. hits") + xlab("No. peptides") +
  theme(legend.position = "none")
ggsave("General_Y.score1/Personalized_MeanMed.svg", width = 5, height = 3)
