# script to load NSBP data and perform metaMDS
# Data is from ICES North Sea Benthos Project
# 
# Author: evberghe
###############################################################################

library(vegan) # for vegdist and metaMDS
library(reshape2) # for the cast functions
curdir <- "/home/evberghe/workspace/r_scratch/vegan"
list.files(curdir)

# read the raw data, and add a column which combines station and replicate names in a single sample name
df_raw <- read.table(paste(curdir, "distribution_a.txt", sep="/"), header=TRUE, sep=",", quote="\"")
df_raw <- cbind(df_raw, sample=paste(df_raw$Station, df_raw$Replicate, sep="_"))

# select the relevant columns in a new dataframe
# 'all' because this contains a row for each of the available data points
# rename the columns tosomething that makes sense outside the context of NSBP
df_all <- df_raw[, c("id", "sample", "Cleanname", "individuals")]
names(df_all) <- c("id", "sample", "taxon", "count")

# Resolve some data issues
# 'presence' for which there was no count is coded as '-1', change to 1
# some counts were reconstructed from density and sample size, and are not integer
# some of the reconstructed counts are between 0 and 0.5, so round would create counts of zero
df_all$count[(df_all$count>0 & df_all$count<=0.5) | df_all$count==-1] <- 1
df_all$count <- round(df_all$count)
df_all <- df_all[df_all$count!=0,]

# subsample
totrecs <- length(df_all$id)
subsamplesize <- 10000
subsample <- sample.int(totrecs, subsamplesize)
df <- df_all[subsample, ]
# 'sample' and 'taxon' are factors; subsampling will have caused some factor levels
# to be present in the factor object, but not in the data
# reapply factor() to drop unused factor levels
df$sample <- factor(df$sample)
df$taxon <- factor(df$taxon)

# vegdist wants a table, not a list; reformat
ar <- acast(df, sample~taxon, sum, value.var="count")
# look at the resulting data table in a pop-up grid
# edit(ar)

#run the stuff
nsbp_dist_bray <- vegdist(ar)
nsbp_mds_bray <- metaMDS(nsbp_dist_bray)

