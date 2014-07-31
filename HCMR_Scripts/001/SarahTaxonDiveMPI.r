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
setwd("/home/patkos/evaluation/001")

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


} #end if (myrank==0)

# TP: Broadcast the two matrices and the step
agg.all <- bcast(agg, rank.source = 0)
add.all <- bcast(add, rank.source = 0)

step <- bcast(mtx.step, rank.source = 0)

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


#comm.print(dim(out),all.rank=T)


###############
##### taxondive
###############
# Here starts the code for taxondive
#resu <- taxondive(comm2, taxdis, match.force = FALSE)

comm <- NULL

# TP: Only the master node performs these tasks
if (myrank == 0) {
comm2<- read.table("matrixSpecies_noHeader_Percent.csv", header = FALSE, sep=",")

#the community file has to have taxon names as columns, and NO categories in the first row, so it has to be transposed. Since it cannot be transposed into a numerical matrix (it contains text as well and results in a character matrix) it has to be done stepwise: 

# transpose only the numerical parts, leaving the first column as it is
comm<-t(comm2[,-1])
#convert it into a data frame
comm<-data.frame(comm)
#assign column names to the new data frame by taking the first column from the original data frame 
names(comm) <- comm2[,1] 

comm <- as.matrix(comm)

    # These are the steps we want to parallelize
    #del <- dstar <- dplus <- Ed <- Edstar <- edplus <- NULL
    #del <- apply(comm, 1, function(x) sum(as.dist(outer(x, x)) * dis))
    #dstar <- apply(comm, 1, function(x) sum(dis * (xx <- as.dist(outer(x, x))))/sum(xx))
    #rs <- rowSums(comm)
    #del <- del/rs/(rs - 1) * 2
    #cs <- colSums(comm)
    #tmp <- sum(as.dist(outer(cs, cs)) * dis)
    #Ed <- tmp/sum(cs)/sum(cs - 1) * 2
    #Edstar <- tmp/sum(cs)/(sum(cs) - 1) * 2
  
}


#### del and dstar
#del <- apply(comm, 1, function(x) sum(as.dist(outer(x, x)) * dis))
del <- NULL
#dstar <- apply(comm, 1, function(x) sum(dis * (xx <- as.dist(outer(x, x))))/sum(xx))
dstar <- NULL

Lambda<-NULL

#1. Broadcast comm to all processors
comm.all <- bcast(comm, rank.source = 0)

for (i in 1:nrow(comm.all)) {
  #2. Calculate partial outer product
  # Each processor has a number of rows of the general outer product NxN matrix
  temp1 <- outer(comm.all[i,start:end], comm.all[i,])

   # Set the diagonal as zero, so that the sums that follow are for dist matrix
   for (j in 1:(end-start+1) )
      temp1[j,(start+j-1)] = 0
  
  # What we need is element by element multiplication, not matrix multiplication! Easy stuff
  sum.dist.outer <- sum(temp1*out)/2

  del[i] <- reduce(sum.dist.outer, op = "sum")
  
  
  # Lambda
  sum.Lambda <- sum(temp1*(out^2))/2
  Lambda[i] <- reduce(sum.Lambda, op = "sum")
  
  # dstar
  xx <- reduce(sum(temp1), op = "sum")/2 

  dstar[i] <- reduce(sum.dist.outer, op = "sum")/xx
}
# CORRECT! Now del and dstar are only on processor 0 and are the same as the serial one
# Notice that we lose some digits though!

#  if (myrank==0) print("del, dstar")
#  comm.print(del, rank.source =0)    
#  comm.print(dstar, rank.source =0)    


#### del and dstar
    rs <- rowSums(comm.all)
    del <- del/rs/(rs - 1) * 2
    cs <- colSums(comm.all)

    # Again, temp1 is distributed among the processors
  temp1 <- outer(cs[start:end], cs)

   # Set the diagonal as zero, so that the sums that follow are for dist matrix
   for (j in 1:(end-start+1) )
      temp1[j,(start+j-1)] = 0
  
  # What we need is element by element multiplication, not matrix multiplication! Easy stuff
  sum.dist.outer <- sum(temp1*out)/2

  tmp <- reduce(sum.dist.outer, op = "sum")
  Ed <- tmp/sum(cs)/sum(cs - 1) * 2
  Edstar <- tmp/sum(cs)/(sum(cs) - 1) * 2
  
# CORRECT! Now Ed and Edstar are only on processor 0 and are the same as the serial one
# Notice that we lose some digits though!

#  if (myrank==0) print("Ed, Edstar")
#  comm.print(Ed, rank.source =0)    
#  comm.print(Edstar, rank.source =0)        

#  comm <- ifelse(comm > 0, 1, 0)
# Probably not needed, its when the matrix is not binary

#  dplus <- apply(comm, 1, function(x) sum(as.dist(outer(x,x)) * dis))
#  Lambda <- apply(comm, 1, function(x) sum(as.dist(outer(x,x)) * dis^2))
#   m <- rowSums(comm)
#  dplus <- dplus/m/(m - 1) * 2
#  Lambda <- Lambda/m/(m - 1) * 2 - dplus^2

dplus <- del
Lambda <- Lambda/rs/(rs - 1) * 2 - dplus^2

#  if (myrank==0) print("dplus, Lambda")
#  comm.print(dplus, rank.source =0)    
#  comm.print(Lambda, rank.source =0)    

#  S <- attr(dis, "Size")
#  omebar <- sum(dis)/S/(S - 1) * 2
#  varome <- sum(dis^2)/S/(S - 1) * 2 - omebar^2
#  omei <- rowSums(as.matrix(dis))/(S - 1)
  S <- ncol(out)

  sum.out <- sum(out)/2
  #omebar <- NULL
  omebar <- reduce(sum.out, op = "sum")
  omebar <- omebar/S/(S - 1) * 2
  
  sum.out <- sum(out^2)/2
  varome <- reduce(sum.out, op = "sum")
  varome <- varome/S/(S - 1) * 2 - omebar^2

#  if (myrank==0) print("omebar, varome")
#  comm.print(omebar, rank.source =0)    
#  comm.print(varome, rank.source =0)    
  
#  omei <- rowSums(as.matrix(dis))/(S - 1)
#  varomebar <- sum(omei^2)/S - omebar^2
#  vardplus <- 2 * (S - m)/(m * (m - 1) * (S - 2) * (S - 3)) * 
#    ((S - m - 1) * varome + 2 * (S - 1) * (m - 2) * varomebar)

  omei.partial <- 0
  for (i in 1:nrow(out))
    omei.partial <- omei.partial + (sum(out[i,])/(S-1))^2

  omei <- reduce(omei.partial, op = "sum")
  varomebar <- omei/S - omebar^2
  vardplus <- 2 * (S - rs)/(rs * (rs - 1) * (S - 2) * (S - 3)) * ((S - rs - 1) * varome + 2 * (S - 1) * (rs - 2) * varomebar)    

#  if (myrank==0) print("varomebar, vardplus")
#  comm.print(varomebar, rank.source =0)    
#  comm.print(vardplus, rank.source =0)    

  Taxondive.out <- list(Species = rs, D = del, Dstar = dstar, Lambda = Lambda, 
              Dplus = dplus, sd.Dplus = sqrt(vardplus), SDplus = rs * 
                dplus, ED = Ed, EDstar = Edstar, EDplus = omebar)
  class(Taxondive.out) <- "taxondive"



# Stop the clock
comm.print("Full program execution time for each processor",rank.source = 0)
comm.print(proc.time() - fullProgramTimer, all.rank=T)

comm.print("Dimension of the biggest matrix needed during calculation and distributed among the processors",rank.source = 0)
comm.print(dim(temp1),all.rank=T)

comm.print("Size of the biggest matrix needed during calculation and distributed among the processors (Mb)",rank.source = 0)
comm.print(object.size(temp1)/1048600,all.rank=T)

  
comm.print("Taxondive values (gathered in processor 0)",rank.source = 0)


  #if (myrank==0) print("Species")
  comm.print(Taxondive.out[1], rank.source =0)    
  #if (myrank==0) print("D")
  comm.print(Taxondive.out[2], rank.source =0)    
  #if (myrank==0) print("Dstar")
  comm.print(Taxondive.out[3], rank.source =0)    
  #if (myrank==0) print("Lambda")
  comm.print(Taxondive.out[4], rank.source =0)    
  #if (myrank==0) print("Dplus")
  comm.print(Taxondive.out[5], rank.source =0)    
  #if (myrank==0) print("sd.Dplus")
  comm.print(Taxondive.out[6], rank.source =0)    
  #if (myrank==0) print("SDplus")
  comm.print(Taxondive.out[7], rank.source =0)    
  #if (myrank==0) print("ED")
  comm.print(Taxondive.out[8], rank.source =0)    
  #if (myrank==0) print("EDstar")
  comm.print(Taxondive.out[9], rank.source =0)
  #if (myrank==0) print("EDplus")
  comm.print(Taxondive.out[10], rank.source =0)    


 finalize()
