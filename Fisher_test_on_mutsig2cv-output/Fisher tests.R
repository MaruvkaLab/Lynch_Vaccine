# Author: Hagay Ladany
# Date: 8.dec.2024
# Load data
Merged.counts <- read.csv("Merged.counts.csv") # Output summary of SNV calls - counts of significant genes found using MutSig2CV https://hub.docker.com/r/tejasrj/mutsig2cv
GENES <- colnames(Merged.counts)[-(1:8)] # List of genes
# ------------------------------------------------------------------------------
# Fisher tests: Lynch Adenoma vs. Carcinoma
# ------------------------------------------------------------------------------
fisher.aden.carc.Lynch <- list()
for (GENE in GENES){
  TEMP <-tapply(Merged.counts[which(Merged.counts$class.=="Lynch"),GENE],Merged.counts[which(Merged.counts$class.=="Lynch"),"Histology"],function(X) X[!is.na(X)])
  TEMP.mat <- matrix(data = c(true.aden=sum(TEMP[["ADENOMA"]]>0),
                            false.aden=sum(TEMP[["ADENOMA"]]==0),
                            true.carc=sum(TEMP[["CARCINOMA"]]>0),
                            false.carc=sum(TEMP[["CARCINOMA"]]==0)),
                   ncol = 2,byrow = F,dimnames = list(c("Mut","WT"),c("Aden.","Carc.")))
  fisher.aden.carc.Lynch[[paste0("less"     ,"_",GENE)]]  <-  fisher.test(TEMP.mat,alternative = "less"     )$p.value
  fisher.aden.carc.Lynch[[paste0("two.sided","_",GENE)]]  <-  fisher.test(TEMP.mat,alternative = "two.sided")$p.value
}
fisher.aden.carc.Lynch <- data.frame(p.val=unlist(fisher.aden.carc.Lynch))
fisher.aden.carc.Lynch <- data.frame(tapply(fisher.aden.carc.Lynch, sub("_.*","",row.names(fisher.aden.carc.Lynch)), rbind),row.names = GENES)
fisher.aden.carc.Lynch$q.val.less       <-  p.adjust(fisher.aden.carc.Lynch$less,method = "BH")
fisher.aden.carc.Lynch$q.val.two.sided  <-  p.adjust(fisher.aden.carc.Lynch$two.sided,method = "BH")
Aden.Carc.fisher <- rownames(fisher.aden.carc.Lynch)[fisher.aden.carc.Lynch$q.val.less<.1]
cat("Fisher tests: Lynch Adenoma vs. Carcinoma: ",Aden.Carc.fisher)
# ------------------------------------------------------------------------------
# Fisher tests: Lynch Carcinoma vs. Sporadic-MSI
# ------------------------------------------------------------------------------
fisher.Lync.carconomas.Spo  <-  list()
for (GENE in GENES){
  TEMP <- tapply(Merged.counts[Merged.counts$Histology=="CARCINOMA",GENE],Merged.counts[Merged.counts$Histology=="CARCINOMA","class."],function(X) X[!is.na(X)])
  TEMP.mat <- matrix(data = c(true.aden=sum(TEMP[["Lynch"]]>0),
                            false.aden=sum(TEMP[["Lynch"]]==0),
                            true.carc=sum(TEMP[["Sporadic"]]>0),
                            false.carc=sum(TEMP[["Sporadic"]]==0)),
                   ncol = 2,byrow = F,dimnames = list(c("Mut","WT"),c("Lynch","Sporadic")))
  fisher.Lync.carconomas.Spo[[paste0("greater"     ,"_",GENE)]] <- fisher.test(TEMP.mat,alternative = "greater")$p.value
  fisher.Lync.carconomas.Spo[[paste0("less"     ,"_",GENE)]] <- fisher.test(TEMP.mat,alternative = "less")$p.value
  fisher.Lync.carconomas.Spo[[paste0("two.sided","_",GENE)]] <- fisher.test(TEMP.mat,alternative = "two.sided")$p.value
}
fisher.Lync.carconomas.Spo <- data.frame(p.val=unlist(fisher.Lync.carconomas.Spo))
fisher.Lync.carconomas.Spo <- data.frame(tapply(fisher.Lync.carconomas.Spo, sub("_.*","",row.names(fisher.Lync.carconomas.Spo)), rbind),row.names = GENES)
fisher.Lync.carconomas.Spo$q.val.greater <-  p.adjust(fisher.Lync.carconomas.Spo$greater,method = "BH")
fisher.Lync.carconomas.Spo$q.val.less <-  p.adjust(fisher.Lync.carconomas.Spo$less,method = "BH")
fisher.Lync.carconomas.Spo$q.val.two.sided <-  p.adjust(fisher.Lync.carconomas.Spo$two.sided,method = "BH")
LynCarc.Sporadic.fisher <- rownames(fisher.Lync.carconomas.Spo)[fisher.Lync.carconomas.Spo$q.val.two.sided<.1]
cat("Fisher tests: Lynch Carcinoma vs. Sporadic-MSI: ",LynCarc.Sporadic.fisher)
