#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(stringi)
library(optparse)


###############################
##       SET VARIABLES       ##


option_list = list(
		make_option(c("-w", "--workPATH"), type="character", default=NULL, 
				help="Main working directory, where list of individual sample (full) names are located.", metavar="character"),
		make_option(c("-i", "--inputPATH"), type="character", default=NULL, 
				help="folder name where to find input files (should be a subfolder in workPATH or its path be provided relative to workPATH)", metavar="character"),
		make_option(c("-o", "--outputPATH"), type="character", default="R_z-link_SCAFF", 
				help="folder name where to export output files (will be created in workPATH). [default= %default]", metavar="character"),
		make_option(c("-s", "--species"), type="character", default=NULL, 
				help="species name", metavar="character")
		)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$workPATH)){
	print_help(opt_parser)
	stop("The main working directory needs to be supplied.\n", call.=FALSE)
} else if (is.null(opt$inputPATH)){
	print_help(opt_parser)
	stop("The input file name needs to be supplied.\n", call.=FALSE)
} else if (is.null(opt$species)){
	print_help(opt_parser)
	stop("The output file name needs to be supplied.\n", call.=FALSE)
}

workDIR = opt$workPATH				# workDIR = "/group/sbs_ssin/stella/albatross/05_reseq"
iDIR = paste0(opt$workPATH,"/",opt$inputPATH)	# iDIR = paste0(workDIR, "/05_GATK")
oDIR = paste0(opt$workPATH,"/",opt$outputPATH)	# oDIR = paste0(workDIR, "/06_ANGSD")
SPECIES = opt$species				# SPECIES = "BFAL" / "LAAL"



###################################
##     SET ENVIRONMENT PATHS     ##

setwd(workDIR)
if( ! dir.exists(oDIR) ) { dir.create(oDIR) }


########################################
###          Start analyses          ###

samples = readLines(paste0(workDIR,"/list_reseq_",SPECIES))
samples = strsplit(samples, "_(?=[^_]+$)", perl=TRUE)
samples = unlist(lapply(samples, `[`, 1))

if(SPECIES=="BFAL") { samples = c(samples, "BFAL_220bp_700") }


# 1 ################

df = data.frame()

s = samples[1]
for(s in samples) {
	infile = read.table(paste0(iDIR,"/idxstats_",s,".txt"))
	infile$sample = paste0(s)
	df = rbind(df,infile)
}

rm(s,infile)

unique(df$sample)
colnames(df) = c("scaffold","length","no_mapped","no_unmapped","sample")


# 2 ################

unmappedread_summary = filter(df, scaffold == "*") # this might be useful later
df = filter(df, scaffold != "*") # removes the last line of each infile (sum of unmapped reads)
tbl_df(df)


# 3 ################


# make new column of coverage estimate
df = mutate(df, cov_est = 125*no_mapped/length) # 125 = read lengths (bp) ??

# make new column of median coverage per sample
   # (the mean was not used because is right-skewed by some large values)
df <- df %>%
  group_by(sample) %>%
  mutate(med_cov = median(cov_est))

# make new column of standardised median coverage per sample
df <- df %>%
  group_by(sample) %>%
  mutate(scaf_cov_stand = cov_est/med_cov)


# 4 ################

#make a new table of median coverage per sample
samp_cov <- df %>%
  group_by(sample) %>%
  summarise(first(med_cov))
write.csv(samp_cov, file=paste0(oDIR,"/", SPECIES,"_per-samp_cov_est.csv"))


# 5 ################


# include the sex designation, acertained by sexing PCR:

if (SPECIES=="BFAL") {

	df$sex <- ifelse(df$sample=="BFAL_1_BFAL676","FEMALE", 
		  ifelse(df$sample=="BFAL_2_3","FEMALE", 
          	  ifelse(df$sample=="BFAL_4_HEW236","FEMALE", 
          	  ifelse(df$sample=="BFAL_10_20","FEMALE", 
          	  ifelse(df$sample=="BFAL_11_1","FEMALE", 
          	  ifelse(df$sample=="BFAL_14_BFAL680","FEMALE", 
          	  ifelse(df$sample=="BFAL_15_18","FEMALE","MALE")))))))
	# Males below:
	#          ifelse(df$sample=="BFAL_12_BFAL674 \
	#          ifelse(df$sample=="BFAL_13_HEW052 \
	#          ifelse(df$sample=="BFAL_3_BFAL698 \
	#          ifelse(df$sample=="BFAL_5_HEW077 \
	#          ifelse(df$sample=="BFAL_6_HEW182 \
	#          ifelse(df$sample=="BFAL_7_BFAL691 \
	#          ifelse(df$sample=="BFAL_8_HEW072 \
	#          ifelse(df$sample=="BFAL_9_4 \
	#          ifelse(df$sample=="BFAL_220bp_700

} else if (SPECIES=="LAAL") {

	df$sex <- ifelse(df$sample=="LAAL_18_HEW169","FEMALE",
          	  ifelse(df$sample=="LAAL_19_HEW173","FEMALE",
          	  ifelse(df$sample=="LAAL_21_HEW179","FEMALE",
          	  ifelse(df$sample=="LAAL_22_HEW188","FEMALE",
          	  ifelse(df$sample=="LAAL_23_HEW129","FEMALE",
          	  ifelse(df$sample=="LAAL_24_HEW133","FEMALE",
          	  ifelse(df$sample=="LAAL_29_HEW218","FEMALE","MALE")))))))
	# Males below:
	#          ifelse(df$sample=="LAAL_16_HEW116 \
	#          ifelse(df$sample=="LAAL_17_HEW162 \
	#          ifelse(df$sample=="LAAL_20_HEW176 \
	#          ifelse(df$sample=="LAAL_25_HEW143 \
	#          ifelse(df$sample=="LAAL_26_HEW158 \
	#          ifelse(df$sample=="LAAL_27_HEW192 \
	#          ifelse(df$sample=="LAAL_28_HEW204 \
	#          ifelse(df$sample=="LAAL_30_HEW231 \
	#          ifelse(df$sample=="LAAL_220bp_HEW162

}

df$sex = as.factor(df$sex) #@@@@ ADDITIONAL LINE FOR t.test() 

df = separate(df, scaffold, c("del","scaf"), sep="_", remove=FALSE) #require
df$del = NULL


# remove scaf 2671, which is the only scaf that has zero reads for some samples (why??)
#df<-filter(df, scaf!=2671) 
#@@@@ THIS ISSUE DISAPPEARED ONCE "df = filter(df, length != 0)" has been replaced by "df = filter(df, scaffold != "*")" 



# 6 ################


# summarise the data by gender for each scaffold
df_new <- df %>%
  group_by(scaffold)  %>%
  summarise(scaf_length = first(length),
            mean_depth_male = mean(scaf_cov_stand[which(sex=="MALE")]),
            med_depth_male = median(scaf_cov_stand[which(sex=="MALE")]),
            var_depth_male = var(scaf_cov_stand[which(sex=="MALE")]),
            n_male = n_distinct(scaf_cov_stand[which(sex=="MALE")]),
            mean_depth_female = mean(scaf_cov_stand[which(sex=="FEMALE")]),
            med_depth_female = median(scaf_cov_stand[which(sex=="FEMALE")]),
            var_depth_female = var(scaf_cov_stand[which(sex=="FEMALE")]),
            n_female = n_distinct(scaf_cov_stand[which(sex=="FEMALE")]))

df_new = filter(df_new, scaffold != "*") #@@@@ WHY "*" APPEARED BACK ???

df_new = separate(df_new,scaffold, c("del","scaf"),sep="_",remove=FALSE)
df_new$del = NULL
df_new$scaf = as.numeric(df_new$scaf)

df_new = arrange(df_new,scaf)
df_new <- df_new %>%
  group_by(scaffold) %>%
  mutate(ratio_mf = mean_depth_male/mean_depth_female)


numb_scafs = length(df_new$scaffold)
#Add a new column in df_new for the p-values, calculated from df 

for(i in c(1:length(df_new$scaffold))){
  df_new$pvalue[i] = t.test(scaf_cov_stand ~ sex, data = df[which(df$scaf==i-1),])$p.value
} 


##?? Something wrong with the p-value? Input data?
write.csv(df_new, file=paste0(oDIR,"/", SPECIES,"_df_new-",numb_scafs,"scafs_coverage.csv"))



# 7 ################

pdf(paste0(oDIR,"/BFAL_plot_scaffold_sizes.pdf"))

  # select a scaffold size cutoff (bp)
  plot(df_new$scaf, df_new$scaf_length)
  #min_scaf = 1e5
  min_scaf = 1e4
  abline(h=min_scaf)

dev.off()


# Remove short scaffolds (based on cutoff above)
#df_size_filt = filter(df_new, scaf_length > min_scaf)
#unique(df_size_filt$scaf[which(df_size_filt$mean_depth_male>3)]) # 6 scafs have > 3 times the standardised coverage. Determined that none of these are likely chr 5 candidates - they have v. low Het, and v. high coverage (repeat regions?)



#df_size_filt = filter(df_size_filt, mean_depth_male < 3)

#numb_scafs = length(df_size_filt$scaffold)
# Add a new column in df_new for the p-values, calculated from df 
#for(i in c(1:numb_scafs)){
#  df_size_filt$pvalue[i] = t.test(scaf_cov_stand ~ sex, data = df[which(df$scaf==i-1),])$p.value
#} 
# expect the warnings if using less than the full number of scafs




# checkpoint
#write.csv(df_size_filt, file="/Users/sinyw/Desktop/Harvard postdoc/projects/fairywren project/re-sequencing analyses/analyses scripts/WSFW/idxResults/perscaf_mf_cov_ratio.csv", col.names=TRUE, row.names = FALSE)
#df_size_filt = read_csv(file="/Users/sinyw/Desktop/Harvard postdoc/projects/fairywren project/re-sequencing analyses/analyses scripts/WSFW/idxResults/perscaf_mf_cov_ratio.csv")


write.csv(df_new, file=paste0(oDIR,"/", SPECIES,"_perscaf_mf_cov_ratio.csv"), col.names=T, row.names=F) #@@@@ ADDED LINE


# 8 ################

pdf(paste0(oDIR,"/", SPECIES,"_plot_hist_perscaf_mf_cov_ratio.pdf"))

  # check out the distribution of the per-scaffold m:f coverage ratio. The little peak around 2 likely reveals the z-linked scafs.
  #hist(df_size_filt$ratio_mf, xlim=c(0,3), breaks=50)
  #plot(density(df_size_filt$ratio_mf, xlim=c(0,3)))
  hist(df_new$ratio_mf, xlim=c(0,3), breaks=50)
  abline(v=c(0.6,1.3), col="red", lty = 2)
  abline(v=c(1.8,2.2), col="blue", lty = 2)
  # How to determine these ablines?

dev.off()


########################################################################

#@@@@ CHANGED ALL "df_size_filt" TO "df_new"

# But what about the scafs with ratio <1, 1<ratio<2 and ratio>2. Are they autosome-linked or z-linked scafs, or something weirdd?? Some "arbitrary" cutoff will be required to designate z-linked chroms!

# Make three different BED files for filtering SNPs examined in various analyses:
# 1 - scafs that are very likely autosomal (0.6<r<1.4, p>0.05)
# 2 - scafs that are very likely z-linked (1.6<r>2.4, p<0.05)
# 3 - a "less conservative" autosome bed file?


# Should not use p-value = 0.05 as cut off
#hist(df_size_filt$pvalue, breaks=100)




###
###  1 - autosomal scaffolds  ###

# without pval cut-off
autosomes = filter(df_new, ratio_mf > 0.6 & ratio_mf < 1.3)
sum(autosomes$scaf_length) # [1] 1147463157
write.table(autosomes, file=paste0(oDIR,"/", SPECIES,"_list_autosomes.txt"), quote=F, row.names=F, sep="\t")

scafs_1 = as.data.frame(cbind(as.character(autosomes$scaffold), 0, autosomes$scaf_length))
write.table(scafs_1, file=paste0(oDIR,"/", SPECIES,"_autosome_scafs_1e4.bed")
            , sep="\t", append=F, quote=F, row.names=F, col.names=F)


# with correction for multiple tests (pval= 0.05/2042)
autosomes_pval = filter(df_new, ratio_mf > 0.6 & ratio_mf < 1.3 & pvalue > 0.000024485798237)
sum(autosomes_pval$scaf_length) # [1] 1121702017
write.table(autosomes_pval, file=paste0(oDIR,"/", SPECIES,"_list_autosomes_pval.txt"), quote=F, row.names=F, sep="\t")

scafs_1p = as.data.frame(cbind(as.character(autosomes_pval$scaffold), 0, autosomes_pval$scaf_length))
write.table(scafs_1p, file=paste0(oDIR,"/", SPECIES,"_autosome_scafs_pval_1e4.bed")
            , sep="\t", append=F, quote=F, row.names=F, col.names=F)



###
###  2 - z-linked scaffolds ###

# without pval cut-off
zchrom = filter(df_new, ratio_mf > 1.5 & ratio_mf < 2.2)
sum(zchrom$scaf_length) # [1] 54540034
write.table(zchrom, file=paste0(oDIR,"/", SPECIES,"_list_zchrom.txt"), quote=F, row.names=F, sep="\t")

scafs_2 = as.data.frame(cbind(as.character(zchrom$scaffold), 0, zchrom$scaf_length))
write.table(scafs_2, file=paste0(oDIR,"/", SPECIES,"_zchrom_scafs_1e4.bed")
            , sep="\t", append=F, quote=F, row.names=F, col.names=F)


# with correction for multiple tests (pval= 0.05/2042)
zchrom_pval = filter(df_new, ratio_mf > 1.5 & ratio_mf < 2.2 & pvalue < 0.000024485798237)
sum(zchrom_pval$scaf_length) # [1] 49649825
write.table(zchrom_pval, file=paste0(oDIR,"/", SPECIES,"_list_zchrom_pval.txt"), quote=F, row.names=F, sep="\t")

scafs_2p = as.data.frame(cbind(as.character(zchrom_pval$scaffold), 0, zchrom_pval$scaf_length))
write.table(scafs_2p, file=paste0(oDIR,"/", SPECIES,"_zchrom_scafs_pval_1e4.bed")
            , sep="\t", append=F, quote=F, row.names=F, col.names=F)




########################################################################


# NB: autosome_scafs_bed seemed to catch many z-linked scafs (min size = 100000bp[?], ratio: 0.7-1.2)
#autos2 <- as.data.frame(cbind(as.character(autosomes$scaffold), 0, autosomes$scaf_length)) 
#write.table(autos2, file=paste0(oDIR,"BFAL_autosome_scafs_1e5.bed"), sep="\t", append=FALSE, quote=FALSE, row.names = FALSE, col.names = c("chrom","chromStart","chromEnd"))












