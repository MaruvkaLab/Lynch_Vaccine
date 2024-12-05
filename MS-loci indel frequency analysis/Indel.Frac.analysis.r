# MMRd Colorectal Neoplasia Statistical Analysis Script
# This script performs pairwise Fisher's Exact Tests for mutation data.
# The script was adjusted for public availability and was altered.
# Author: Hagay Ladany
# Date: 4.Dec.2024
# ================================
# Load Required Libraries
# ================================
library(data.table)
library(readxl)
library(biomaRt)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
# ================================
# Define Helper Functions
# ================================
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# ================================
# Load Data
# ================================
MetaData <- read_excel("Supplementary Table 1 - Metadata.xlsx", na = "NA")  # Read the Excel metadata file
load("Relevant_Pos.snvMatched.mini.Rdata") #load loci data
# ================================
# Function to get gene symbol given chromosome and start position
# ================================
ensembl <- useMart(host = "https://grch37.ensembl.org","ensembl", dataset = "hsapiens_gene_ensembl")
get_gene_symbol <- function(chr, start) {
  # Perform query
  result <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position','gene_biotype'),
                  filters = c('chromosome_name', 'start', 'end'),
                  values = list(chr, start, start),
                  mart = ensembl)
  result <- result[result$hgnc_symbol!="",]
  
  # Check if result is empty
  if (nrow(result) == 0) {
    return("NA")
  } else {
    return(result[result$gene_biotype=="protein_coding","hgnc_symbol"][1])
  }
}
# ================================
# Function generating a volcano plot and statistics.
# ================================
create_scatter <- function(DATASET,SORT_BY,CAT1,CAT2,PLOT_TT="",PLOT_CAPS=SNV.MATCHING_STAT,ADD_GENE_NAMES=F) {
  cat1 <- firstup(tolower(CAT1))
  cat2 <- firstup(tolower(CAT2))
  DATASET      <- DATASET[DATASET[,SORT_BY]%in%c(CAT1,CAT2),]
  sample_t     <- table(MetaData[MetaData$`Paper ID`%in%unique(DATASET$paper_id),SORT_BY])
  TEMP         <- as.data.frame(table(DATASET[,c("IDs","CALLed",SORT_BY)]))
  if(length(CAT2)>1){cat2<-paste0(CAT2,collapse = "-");
                     TEMP[,SORT_BY]<-as.character(TEMP[,SORT_BY]);
                     TEMP[,SORT_BY][TEMP[,SORT_BY]!=CAT1] <- cat2;
                     TEMP[,SORT_BY]<-as.factor(TEMP[,SORT_BY])
  }
  TEMP2 <- reshape2::dcast(TEMP, IDs~...,value.var = "Freq",fun.aggregate=sum)
  TEMP3 <- cbind(TEMP2[,c(2,4)]/rowSums(TEMP2[,c(2,4)]),TEMP2[,c(3,5)]/rowSums(TEMP2[,c(3,5)]))
  colnames(TEMP3) <- paste0(colnames(TEMP3),"_frac")
  TEMP2 <- data.frame(TEMP2,SUM=rowSums(TEMP2[,-1]),TEMP3)
  colnames(TEMP2) <- sub("FALSE", "F",colnames(TEMP2))
  colnames(TEMP2) <- sub("TRUE", "T",colnames(TEMP2))
  MAT <- data.frame(p.val=tapply(TEMP2[,2:5],TEMP2$IDs,function (X) matrix(unlist(X),nrow = 2,dimnames = list(c("Un-called","Called"),c(cat1,cat2)),byrow = T)))
  Extream <- lapply(MAT[[1]],function(X) matrix(c(0,colSums(X),0),ncol=2,dimnames = list(c("Un-called","Called"),c("Adenoma","Carcinoma"))))
  Extream.rez <- data.frame(p.val=sapply(Extream, function (X) fisher.test(X)$p.value))
  Extream.rez$q.val <- p.adjust(Extream.rez$p.val,method = "BH")>.1
  TEMP.rez <- data.frame(p.val=sapply(MAT[[1]], function (X) fisher.test(matrix(X,nrow = 2))$p.value))
  rez.comparison <- cbind(TEMP2,TEMP.rez)
  # Remove unwanted tests
  rez.comparison <- rez.comparison[!Extream.rez$q.val,]
  rez.comparison$q.val <- p.adjust(rez.comparison$p.val,method = "BH")
  rez.comparison$PASS <- rez.comparison$q.val< 0.1
  if(!all(!rez.comparison$PASS)){temp <- table(rez.comparison$PASS)/sum(table(rez.comparison$PASS))}else{temp <- 0 }
  if(ADD_GENE_NAMES){
    # Order by significant
    rez.comparison  <-  rez.comparison[order(rez.comparison$p.val,decreasing = F),]
    # Get gene names 
    TempGN<-data.frame(do.call(rbind,strsplit(head(as.character(rez.comparison[1:10,"IDs"]),n=10),"_"))[,1:2])
    rez.comparison$GENE_name<-NA
    rez.comparison$GENE_name[1:10] <- apply(TempGN,1,function(X) get_gene_symbol(chr = X[["X1"]],start = as.numeric(X[["X2"]])))
  }
  # Create plot
  Graph_volcan<-ggplot(data = NULL,aes(x=log10(rez.comparison[,paste0("T_",CAT1,"_frac")]/rez.comparison[,paste0("T_",paste0(CAT2,collapse = "."),"_frac")]),
                                       y=-log10(rez.comparison$p.val)))+
    geom_point(size=.25*pracma::nthroot(rez.comparison$SUM,n = 2),alpha=.33,shape=21,fill=ifelse(rez.comparison$PASS,"red3","lightblue"))+
    labs(title = PLOT_TT,
         subtitle = paste0(cat2,"  n=",sum(sample_t[names(sample_t)!=CAT1]),";  ",cat1, if(SORT_BY=="Lynch_germLine_mut"){toupper(cat1)},"  n=",sample_t[[CAT1]]),
         caption = PLOT_CAPS)+
    ylab(expression(-log[10](p.value)))+
    xlab(bquote("MS-Indel Fractions log"[10] * "(" * .(if(SORT_BY=="Lynch_germLine_mut"){toupper(cat1)}else{cat1}) * "/" * .(cat2) * ")"))+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    geom_text_repel(aes(label=rez.comparison$GENE_name),size=2)

  return(list(Raw.table=rez.comparison,Graph_volcan=Graph_volcan))
}

#===================================
# View results:
#===================================
#Lynch Adenoma Vs Carcinoma (without MSH6) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subset samples for analysis
Tmp.L <- Relevant_Pos.snvMatched.mini[-which(Relevant_Pos.snvMatched.mini$Lynch_germLine_mut=="MSH6"),]
Tmp.L <- Tmp.L[which(Tmp.L$class=="Lynch"),]
Adenoma.Carcinoma.noMSH6      <- create_scatter(DATASET = Tmp.L,
                                                SORT_BY = "Histology",
                                                CAT1 = "ADENOMA",
                                                CAT2 = "CARCINOMA",
                                                PLOT_TT = "Lynch-Adenoma vs Lynch-Carcinomas\n(without MSH6 samples)",
                                                PLOT_CAPS = "",
                                                ADD_GENE_NAMES = T)
# view volcano
Adenoma.Carcinoma.noMSH6$Graph_volcan
# view test results
Adenoma.Carcinoma.noMSH6$Raw.table

#Lynch Carcinoma Vs Sporadic-MSI (without MSH6) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subset samples for analysis
Tmp <- Relevant_Pos.snvMatched.mini[which(Relevant_Pos.snvMatched.mini$Histology=="CARCINOMA"),]
Tmp <- Tmp[-which(Tmp$Lynch_germLine_mut=="MSH6"),]
Sporadic_Vs_Lynch_Carcinoma.noMSH6 <- create_scatter(
  DATASET = Tmp,
  SORT_BY = "class.",
  CAT1 = "Lynch",
  CAT2 = "Sporadic",
  PLOT_TT = "Lynch-Carinomas Vs Sporadic-MSI\n(without MSH6 samples)",
  PLOT_CAPS = "",
  ADD_GENE_NAMES = T)
# view volcano
Sporadic_Vs_Lynch_Carcinoma.noMSH6$Graph_volcan
# view test results
Sporadic_Vs_Lynch_Carcinoma.noMSH6$Raw.table
