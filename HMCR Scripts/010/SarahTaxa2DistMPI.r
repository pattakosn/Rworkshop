#Script to calculate Taxonomic Distinctness
# An attempt to parallelize Sarah's code with pbdMPI.

#clean previous results
rm(list=ls())


# Start the clock
fullProgramTimer <- proc.time()

# TP: Include library for MPI communication
library(pbdMPI, quiet = TRUE)

init()

# TP: Get the size of the cluster
processors <- comm.size()
myrank <- comm.rank()


# TP: The step will be the number of rows/columns assigned to each node
# TP: The rest are matrices that will be sent to each node
mtx.step <- NULL
agg <- NULL
add <- NULL


#load vegan library
library(vegan)

# TP: Only the master node performs these tasks
if (myrank == 0) {
# set path to files 
setwd("/home/patkos/evaluation/010")

#aggs contains the aggregation file
agg<- read.table("aggSpecies_Percent.csv", header = TRUE, sep=",")


# taxa2dist is the function we want to parallelize
#taxdis <- taxa2dist(agg, varstep=TRUE)
# Lets delve into its code:
#function (x, varstep = FALSE, check = TRUE, labels) 
varstep = TRUE
check = TRUE

  rich <- apply(agg, 2, function(taxa) length(unique(taxa)))
  S <- nrow(agg)
  if (check) {
    keep <- rich < S & rich > 1
    rich <- rich[keep]
    agg <- agg[, keep]
  }
  i <- rev(order(rich))
  agg <- agg[, i]
  rich <- rich[i]
  if (varstep) {
    add <- -diff(c(nrow(agg), rich, 1))
    add <- add/c(S, rich)
    add <- add/sum(add) * 100
  }
 # else {
 #   add <- rep(100/(ncol(agg) + check), ncol(agg) + check)
 # }
  if (!is.null(names(add))) 
    names(add) <- c("Base", names(add)[-length(add)])
  if (!check) 
    add <- c(0, add)
 
  #########################################
  # TP: Pretty huge matrix... The big dataset generates 18.5 billion elements (168931)
  #########################################
   addTemp <- add
  #out <- matrix(add[1], nrow(agg), nrow(agg))

  
  
# TP: Only matrices can be sent through pbdMPI so we transform 'add' to matrix
add <- matrix(add)
mtx.step <- matrix(ceiling(nrow(agg)/processors),1,1)


}

# TP: Broadcast the two matrices and the step
agg.all <- bcast(agg, rank.source = 0)
add.all <- bcast(add, rank.source = 0)

step <- bcast(mtx.step, rank.source = 0)


# Start the clock
mainCalcTimer <- proc.time()

# TP: Now, each node has a different set of data to work with
start <- myrank*step+1
end <- start+step-1
# TP: The last node has less elements that the rest, i.e., the remaining elements
if (myrank == processors-1) {end<-nrow(agg.all)}

#comm.print(object.size(agg.all), all.rank = TRUE)
#comm.print(start, all.rank = TRUE)
#comm.print(end, all.rank = TRUE)

#out <- matrix(add.all[1], nrow(agg.all), (end-start+1))
out <- matrix(add.all[1], (end-start+1), nrow(agg.all))
#comm.print(dim(out), all.rank = TRUE)


# TP: This is the main loop that is performed on each processor separately
for (i in 1:ncol(agg.all)) {
out <- out + add.all[i + 1] * outer(agg.all[start:end, i], agg.all[, i], "!=")
}

#comm.print(out[1,1], all.rank = TRUE)

# TP: Set the diagonal as 0 (when a distance matrix is transformed into matrix, the diagonal is set to zero)
for (j in 1:nrow(out)) {
   out[j,step*myrank+j] <- 0
}

# Stop the clock
comm.print("Time for each processor to perform the main calculation loop",rank.source = 0)
comm.print(proc.time() - mainCalcTimer, all.rank=T)


# Stop the clock
comm.print("Full program execution time for each processor",rank.source = 0)
comm.print(proc.time() - fullProgramTimer, all.rank=T)

comm.print("Dimension of the matrix stored in each processor",rank.source = 0)
comm.print(dim(out),all.rank=T)

comm.print("Size of the matrix stored in each processor (Mb)",rank.source = 0)
comm.print(object.size(out)/1048600,all.rank=T)


####### Check that the result is the same as the serial
####### Remove this code when the data is too big to fit in processor 0
#if (myrank == 0) {
#  outSerial <- matrix(addTemp[1], nrow(agg), nrow(agg))
  
#  for (i in 1:ncol(agg)) {
#      outSerial <- outSerial + addTemp[i + 1] * outer(agg[, i], agg[, i], "!=")
#  }
  
#outSerial <- as.dist(outSerial)
#outSerial <- as.matrix(outSerial)
#}

# Gather everything at processor 0
#allData <- allgather(out)
#outComplete <- allData[[1]]
#for (j in 2:length(allData)) {
#   outComplete <- rbind(outComplete,allData[[j]])
#}
##comm.print(dim(outComplete))
##comm.print(dim(outSerial))
#comm.print(object.size(outComplete))
#comm.print(object.size(outSerial))
##comm.print(outComplete[1,1])
##comm.print(outSerial[1,1])

## We need to make this cyclic transformation some info is attached when made a distance matrix
## We could compare outComplete with outSerial before becoming a dist, but the latter has a value at the diagonal that is removed when it is made a distance
#outComplete <- as.dist(outComplete)
#outComplete <- as.matrix(outComplete)
#comm.print("Are the serial and distributed matrices equal?",rank.source = 0)
#comm.print(all.equal(outComplete,outSerial))


finalize()
