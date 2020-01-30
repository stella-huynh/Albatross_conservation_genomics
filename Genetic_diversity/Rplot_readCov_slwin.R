#!/usr/bin/env Rscript

library(ggplot2)
library(gtools)

setwd("C:/Users/stell/Desktop/Postdoc_HK/01_prj_albatross/05_reseq/06_ANGSD/common/readCov_slwin/new")


SP = $1

infile1 = paste0(SP,"_m_allscafs.posU_matchF.txt")
infile2 = paste0(SP,"_f_allscafs.posU_matchF.txt")

df1 = read.table(infile1,sep="\t")
colnames(df1) = c("REF","Rpos","QRY","totDepth","meanDepth")
df1$value = df1$meanDepth / 8
df1$SEX = "male"

df2 = read.table(infile2,sep="\t")
colnames(df2) = c("REF","Rpos","QRY","totDepth","meanDepth")
df2$value = df2$meanDepth / 7
df2$SEX = "female"


df = rbind(df1,df2)
df = df[which(df$REF!="GRCg6a_MT"),]
df$REF = as.character(df$REF) ; df$REF = as.factor(df$REF)
levels(df$REF) = unlist(lapply(strsplit(levels(df$REF),"_"),`[[`,2))
idx <- grepl("^[^0-9]", levels(df$REF))
df$REF <- factor(df$REF, c(sort(levels(df$REF)[idx]),mixedsort(levels(df$REF)[!idx])))

my_breaks <- function(x) { if (max(x) < 6000) seq(0, 5000, 1000) else seq(0, 15000, 5000) }


pdf(paste0(SP,"_readCov.pdf"),width=24)

ggplot(df,aes(x=Rpos, y=value, col=SEX)) + geom_line() + facet_grid(.~REF, scales="free_x", space="free_x") +
  scale_x_continuous(expand=c(0,0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background=element_rect(fill="white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.spacing = unit(0, "lines"), panel.background = element_rect(color = "grey")) 

dev.off()



d = cbind(df1,df2)
d = d[,c(1:6,11:13)]
colnames(d) = c("REF","Rpos","QRY","totDepth_M","meanDepth_M","value_M","totDepth_F","meanDepth_F","value_F")
d$REF = as.character(d$REF) ; d$REF = as.factor(d$REF)
levels(d$REF) = unlist(lapply(strsplit(levels(d$REF),"_"),`[[`,2))
idx <- grepl("^[^0-9]", levels(d$REF))
d$REF <- factor(d$REF, c(sort(levels(d$REF)[idx]),mixedsort(levels(d$REF)[!idx])))
d$value = d$value_M / d$value_F


pdf(paste0(SP,"_readCov_ratioMF.pdf"),width=24)

ggplot(d,aes(x=Rpos,y=value)) + geom_line() + facet_grid(.~REF, scales="free_x", space="free_x") +
  scale_x_continuous(expand=c(0,0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background=element_rect(fill="white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.spacing = unit(0, "lines"), panel.background = element_rect(color = "grey")) 

dev.off()





