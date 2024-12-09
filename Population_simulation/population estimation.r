############################################################
# Author: Hagya Ladany
# Date: 9-Dec-2024
# Note: This script was modified for publication purposes.
############################################################
library(dplyr)
library(ggridges)
library(ggplot2)
# load population simulation - these were form using the HLA population frequencies of difference races 
# taken from: https://doi.org/10.1038/nbt.3344
# each sample was assign with a set of mutation based on the total frequency and mutation counts distribution found in our cohort.
load("Pseudo_binders.White1.Rdata")
load("Pseudo_binders.Black1.Rdata")
load("Pseudo_binders.Asian1.Rdata")
# merge the data
RACE <- list(Pseudo_bindersW,Pseudo_bindersB,Pseudo_bindersA)
names(RACE) <- c("White","Black","Asian")
# name the samples
for(R in names(RACE)){for(SAMPLE in 1:length(RACE[[R]])){RACE[[R]][[SAMPLE]]$SAMPLE <-SAMPLE}}
# set result object, target administered neo-peptides and binding thresold (70nM)
Rez<-list()
MAXPEP <- 200
BINDING_THRESH <- 70
# estimate simulation binders
for (Rtype in names(RACE)){
# get one binder per mutation
  Temp <- do.call(rbind,RACE[[Rtype]])
  Temp <- Temp[Temp$Binding_score<=BINDING_THRESH,]
  frequency_table <- Temp %>%
    group_by(`MT Epitope Seq`) %>%
    summarise(
      count = n(),
      elements = list(HGVSp),
      samples  = list(SAMPLE))
  
  frequency_table$samples <- lapply(frequency_table$samples,function(X) unique(X))
  frequency_table$nSIZE <- sapply(frequency_table$samples,length)
  frequency_table$elements <- lapply(frequency_table$elements,function(X) unique(X))
  frequency_table$count <- sapply(frequency_table$elements,length)
  frequency_table <- frequency_table[order(decreasing = T,frequency_table$nSIZE),]
  frequency_table <- frequency_table[!duplicated(frequency_table$elements),]
  frequency_table$uniq_Pep <- sapply(frequency_table$elements,function(X) unique(sub(":.*","",X)))
  frequency_table$lenUniPep <- sapply(frequency_table$uniq_Pep,length)
  gc()
  
  Rez[[paste0(Rtype,".hist.data")]]<- data.frame(table(dnn = "Sample",unlist(frequency_table$samples[1:100])),Race=Rtype)
}
# merge results
Binder_race_Hist <- do.call(rbind,Rez)


ggplot(Binder_race_Hist,aes(x=Freq,y=Race,fill=Race,col=Race))+
  geom_density_ridges(rel_min_height=.01)+
  theme_classic()+
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0))+
  scale_x_continuous(limits = c(0,40), expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2,col="purple3") +
  geom_text(data=Binder_race_Hist %>% group_by(Race) %>% 
              summarise(Freq=median(Freq)),
            aes(label= round(Freq)), 
            position=position_nudge(y=.8,x=-2.1), colour="black", size=4.5)+
  # labs(caption = paste0("No. binder < ",BINDING_THRESH,"nM , given the top 100 frequent epitopes in each race\n*Simulated population: n=1000 patients per-race"))+
  scale_fill_manual(values = c("gold2","orange3","aquamarine3"))+
  ylab(NULL)
