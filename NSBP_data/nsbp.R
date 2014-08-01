# script to load NSBP data and perform metaMDS
# Data is from ICES North Sea Benthos Project
# 
# Author: evberghe
###############################################################################

library(vegan)
library(reshape2)
curdir <- "/home/evberghe/workspace/r_scratch/vegan"
list.files(curdir)

df_raw <- read.table(paste(curdir, "distribution_a.txt", sep="/"), header=TRUE, sep=",", quote="\"")
df_raw <- df_raw[, c("id", "Station", "Replicate", "Cleanname", "individuals")]
df_raw <- cbind(df_raw, sample=paste(df_raw$Station, df_raw$Replicate))
df_all <- df_raw[, c("id", "sample", "Cleanname", "individuals")]
names(df_all) <- c("id", "sample", "taxon", "count")
df_all$count[df_all$count==-1] <- 1
df_all$count <- round(df_all$count)

totrecs <- length(df_all$id)
subsamplesize <- 10000
subsample <- sample.int(totrecs, subsamplesize)
df <- df_all[subsample, ]
df$sample <- factor(df$sample)
df$taxon <- factor(df$taxon)
ar <- acast(df, sample~taxon, sum, value.var="count")
edit(ar)

nsbp_dist_bray <- vegdist(ar)
nsbp_mds_bray <- metaMDS(nsbp_dist_bray)

