# Project: STARR - Lynch vaccine 2023
# Script by: Hagay Ladany
# Group: Maruvka Lab, Technion, Israel
# Description: Plotting RNA expression similarity
# Input: RNA counts of STARR samples
# Output: PCA scatter plot
# Libraries
library(data.table)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(EnhancedVolcano)
library(ggpubr)
library(limma)
library(edgeR)
library(clusterProfiler)
library(biomaRt)
# functions
ebsemble_to_hugo <- function(ensemble){
  hugo <- biomaRt::getBM(
    attributes = c("ensembl_gene_id","hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensemble,
    mart = mart)
  return(hugo)
}
# data.preparation ------------------------------------------------------------------------------------------
pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# load col/count data
load("RNA_analysis.Rdata")
#--- limma analysis ------------------------------------------------------------------------------------------
COLDATA <- colDATA[, c("Histopathology", "PAT")]
d0 <- DGEList(rezExp_COUNTS)
# Calculate normalization factors to scale the raw library sizes
d0 <- calcNormFactors(d0)
# Define a cutoff for filtering out low-expressed genes
cutoff<-2
# Identify genes with a maximum counts-per-million (CPM) below the cutoff
drop <- which(apply(cpm(d0), 1, max) < cutoff)
# Filter out the low-expressed genes
d <- d0[-drop,]
# 1.2 Create a design matrix and define contrasts -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
# Create a design matrix for the linear model with the variables of interest
design <- model.matrix(~ 0 + Histopathology + PAT, data = colDATA[, c("Histopathology", "PAT")])
# Define contrasts for the comparisons
contrast.matrix <- makeContrasts(
  normal_vs_carcinoma = HistopathologyCarcinoma - HistopathologyNormal,
  normal_vs_adenoma =   HistopathologyAdenoma   - HistopathologyNormal,
  adenoma_vs_carcinoma =HistopathologyCarcinoma - HistopathologyAdenoma,
  levels = design)
# Define colors for different cultivars
col <- c("orange", "red3", "aquamarine3")[as.numeric(colDATA$Histopathology)]
# 1.3 Voom Transformation -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
# voom transformation -  preparing it for linear modeling    # This includes generating a mean-variance trend plot
d.v <- voom(d, design, plot = F)
#remove batch
d.v.nB <- removeBatchEffect(d.v$E, batch = as.character(colDATA$PAT))
# visualize normalized expression values - 
# 1.4 Fit the Linear Model and Apply Contrasts -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
fit <- lmFit(d.v$E, design)# Fit the linear model to the voom-transformed data
# 1.5 Fit Linear Model and Apply Contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# 1.6 Extracting and Visualizing Results
# Extract the results
results_normal_vs_carcinoma  <- topTable(fit2, coef = "normal_vs_carcinoma" , adjust.method = "fdr", number = Inf)
results_normal_vs_adenoma    <- topTable(fit2, coef = "normal_vs_adenoma"   , adjust.method = "fdr", number = Inf)
results_adenoma_vs_carcinoma <- topTable(fit2, coef = "adenoma_vs_carcinoma", adjust.method = "fdr", number = Inf)
# draw volcano plot------------------------------------------------------------------------------------------------------------------------------------------------------
EnhancedVolcano(results_normal_vs_adenoma, lab = rownames(results_normal_vs_adenoma), x = 'logFC', y = 'P.Value', FCcutoff = 1.5,labSize = 3,subtitle = NULL,legendPosition = "none",caption = NULL,title = "Normal vs Adenoma")+theme(plot.title = element_text(hjust = 0.5))
EnhancedVolcano(results_adenoma_vs_carcinoma, lab = rownames(results_adenoma_vs_carcinoma), x = 'logFC', y = 'P.Value', FCcutoff = 1.5,labSize =3,subtitle = NULL,legendPosition = "none",caption = NULL,title = "Adenoma vs Carcinoma")+theme(plot.title = element_text(hjust = 0.5))
EnhancedVolcano(results_normal_vs_carcinoma, lab = rownames(results_normal_vs_carcinoma), x = 'logFC', y = 'P.Value', FCcutoff = 1.5,labSize = 3,subtitle = NULL,legendPosition = "none",caption = NULL,title = "Normal vs Carcinoma")+theme(plot.title = element_text(hjust = 0.5))
#PCA ---------------------------------------------------------------
pca_result <- prcomp(t(d.v$E), scale. = T) # pca plot
pca_loadings <- pca_result$rotation #Extract PCA Results
pca_scores <- pca_result$x # PCA scores for samples
pca_df <- as.data.frame(pca_scores) # Convert PCA scores to a data frame
pca_df$Sample <- rownames(pca_df)
pca_df$Condition <- colDATA$condition
pca_df$PAT <- colDATA$PAT # Create PCA plot
pca_df$Histopathology <- colDATA$Histopathology
# PCA plot - without batch correction - (figures for Extended data)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Histopathology)) +
  geom_point(size = 3) +
  labs(title = 'PCA of Expression Data',
       x = paste('PC1 (', round(summary(pca_result)$importance[2,1] * 100, 2), '%)', sep = ''),
       y = paste('PC2 (', round(summary(pca_result)$importance[2,2] * 100, 2), '%)', sep = '')) +
  theme_minimal() +
  theme(legend.position = 'right')+
  geom_line(aes(group=PAT),size=.2, alpha=0.5,col="black",linetype=3)+
  geom_text_repel(label=as.numeric(colDATA$PAT),nudge_y = .5,size=2.5)+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=guide_legend(ncol=3,byrow=FALSE,
                             keywidth=2,
                             keyheight=0.1,
                             default.unit="cm"))+
  theme(legend.text = element_text(margin = margin(0, -10, 0, -21)))+
  scale_color_manual(values = c("orangered2", "yellow3", "darkslateblue"))
#PC comparison
#pc1
ggplot(pca_df,aes(x=Histopathology,y=PC1,fill=Histopathology))+
  geom_point()+
  geom_violin(alpha=.2)+
  stat_compare_means(comparisons = list(c("Normal","Adenoma"),c("Normal","Carcinoma"),c("Adenoma","Carcinoma")),method = "wilcox.test",label = "p.signif")+
  theme_bw()+
  theme(legend.position = "none",legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("orangered2", "yellow3", "darkslateblue"))+
  ylab("PC1 value")
#pc2
ggplot(pca_df,aes(x=Histopathology,y=PC2,fill=Histopathology))+
  geom_point()+
  geom_violin(alpha=.2)+
  stat_compare_means(comparisons = list(c("Normal","Adenoma"),c("Normal","Carcinoma"),c("Adenoma","Carcinoma")),method = "wilcox.test",label = "p.signif")+
  theme_bw()+
  theme(legend.position = "none",legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("orangered2", "yellow3", "darkslateblue"))+
  ylab("PC2 value")

# analysis with patient batch corrected (main text)
pca_result <- prcomp(t(d.v.nB), scale. = T) # pca plot
pca_loadings <- pca_result$rotation #Extract PCA Results
pca_scores <- pca_result$x # PCA scores for samples
pca_df <- as.data.frame(pca_scores) # Convert PCA scores to a data frame
pca_df$Sample <- rownames(pca_df)
pca_df$Condition <- colDATA$condition
pca_df$PAT <- colDATA$PAT # Create PCA plot
pca_df$Histopathology <- colDATA$Histopathology
# PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Histopathology)) +
  geom_point(size = 3) +
  labs(title = 'PCA of Expression Data',
       x = paste('PC1 (', round(summary(pca_result)$importance[2,1] * 100, 2), '%)', sep = ''),
       y = paste('PC2 (', round(summary(pca_result)$importance[2,2] * 100, 2), '%)', sep = '')) +
  theme_minimal() +
  theme(legend.position = 'right')+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=guide_legend(ncol=3,byrow=FALSE,
                             keywidth=2,
                             keyheight=0.1,
                             default.unit="cm"))+
  theme(legend.text = element_text(margin = margin(0, -10, 0, -21)))+
  scale_color_manual(values = c("orangered2", "yellow3", "darkslateblue"))
#PC comparison
#pc1
ggplot(pca_df,aes(x=Histopathology,y=PC1,fill=Histopathology))+
  geom_point(size=.5)+
  geom_violin(alpha=.2)+
  stat_compare_means(comparisons = list(c("Normal","Adenoma"),c("Normal","Carcinoma"),c("Adenoma","Carcinoma")),method = "wilcox.test",label = "p.signif")+
  theme_bw()+
  theme(legend.position = "none",legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("orangered2", "yellow3", "darkslateblue"))+
  ylab("PC1 value")
#pc2
ggplot(pca_df,aes(x=Histopathology,y=PC2,fill=Histopathology))+
  geom_point(size=.5)+
  geom_violin(alpha=.2)+
  stat_compare_means(comparisons = list(c("Normal","Adenoma"),c("Normal","Carcinoma"),c("Adenoma","Carcinoma")),method = "wilcox.test",label = "p.signif")+
  theme_bw()+
  theme(legend.position = "none",legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("orangered2", "yellow3", "darkslateblue"))+
  ylab("PC2 value")

# gene set enrichment analysis (GSEA) ----------------------------------------------------
#get hugo symbols
mart <- biomaRt::useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl",host = "https://www.ensembl.org")
TEMP <- ebsemble_to_hugo(sub("\\..*","",rownames(results_normal_vs_carcinoma)))
results_normal_vs_carcinoma$hugo_symbol <- TEMP[match(sub("\\..*","",rownames(results_normal_vs_carcinoma)),TEMP$ensembl_gene_id),"hgnc_symbol"]
TEMP <- ebsemble_to_hugo(sub("\\..*","",rownames(results_normal_vs_adenoma)))
results_normal_vs_adenoma$hugo_symbol <- TEMP[match(sub("\\..*","",rownames(results_normal_vs_adenoma)),TEMP$ensembl_gene_id),"hgnc_symbol"]
TEMP <- ebsemble_to_hugo(sub("\\..*","",rownames(results_adenoma_vs_carcinoma)))
results_adenoma_vs_carcinoma$hugo_symbol <- TEMP[match(sub("\\..*","",rownames(results_adenoma_vs_carcinoma)),TEMP$ensembl_gene_id),"hgnc_symbol"]
# Prepare the input data for GSEA # Select the results of interest, e.g., normal_vs_carcinoma
# Normal vs Carcinoma
results <- results_normal_vs_carcinoma[which(results_normal_vs_carcinoma$hugo_symbol!=""),]
ranked_genes <- results$logFC
names(ranked_genes) <- results$hugo_symbol
ranked_genes <- sort(ranked_genes, decreasing = T)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
N_v_C.H.gseaRes <- fgsea(pathways.hallmark, ranked_genes)
ggplot(N_v_C.H.gseaRes[N_v_C.H.gseaRes$padj<.01,], aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < 0.01)) +  coord_flip() +  labs(x = "Pathway", y = "Normalized Enrichment Score") +  theme_minimal()+ggtitle(NULL)+theme(axis.text.y = element_text(size = 10))
# Normal vs Adenoma
results <- results_normal_vs_adenoma[which(results_normal_vs_adenoma$hugo_symbol!=""),]
ranked_genes <- results$logFC
names(ranked_genes) <- results$hugo_symbol
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
N_v_A.H.gseaRes <- fgsea(pathways = pathways.hallmark, stats = ranked_genes)
ggplot(N_v_A.H.gseaRes[N_v_A.H.gseaRes$padj<.01,], aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < 0.01)) +  coord_flip() +  labs(x = "Pathway", y = "Normalized Enrichment Score") +  theme_minimal()+ggtitle(NULL)+theme(axis.text.y = element_text(size = 10))
# Adenoma vs Carcinoma
results <- results_adenoma_vs_carcinoma[which(results_adenoma_vs_carcinoma$hugo_symbol!=""),]
ranked_genes <- results$logFC
names(ranked_genes) <- results$hugo_symbol
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
A_v_C.H.gseaRes <- fgsea(pathways = pathways.hallmark, stats = ranked_genes)
ggplot(A_v_C.H.gseaRes[A_v_C.H.gseaRes$padj<0.01,], aes(reorder(pathway, NES), NES)) + geom_col(aes(fill = padj < 0.01)) +  coord_flip() +  labs(x = "Pathway", y = "Normalized Enrichment Score") +  theme_minimal() +ggtitle(NULL)+theme(axis.text.y = element_text(size = 10))
