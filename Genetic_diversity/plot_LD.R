#!/usr/bin/Rscript Rscript

library(ggplot2)

FILE=$1

df = read.table(FILE, sep="\t")
colnames(df) = c("sp","chr","minBIN","maxBIN","midLEN","r2")

df$midLEN = (df$minBIN + 500) / 1000000

ggplot(df, aes(x=midLEN, y=r2)) + geom_point(aes(col=sp)) + geom_smooth(aes(col=sp), method="lm") + scale_x_discrete(name="Physical distance (Mb)", limits=seq(0,10,0.5))
  
  
