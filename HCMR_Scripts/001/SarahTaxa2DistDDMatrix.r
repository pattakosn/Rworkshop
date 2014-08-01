#Script to calculate Taxonomic Distinctness

#clean previous results
rm(list=ls())

# Start the clock
fullProgramTimer <- proc.time()

#load vegan library
library(vegan)


# Initialize process grid
library(pbdDMAT, quiet=T)

init.grid()

# set path to files 
setwd("/home/patkos/evaluation/001")


#  if (comm.rank()==0)
#{
# The initialization of data does not need to occur on all nodes. But then we have
# to distribute matrix x, which is inside the outer product! problem..

#aggs contains the aggregation file
agg<- read.table("aggSpecies_Percent.csv", header = TRUE, sep=",")



#comm contains the presence-absence-matrix, but without the header (column names)!
# the delimiter HAS to be ",", otherwise R complains about too long column names (more than 256 bytes, though they are in fact much shorter)
#comm<- read.table("matrixSpecies_noHeader_Percent.csv", header = FALSE, sep=",")


#if needed transform absent variables to 0
#comm[is.na(comm)]<-0


#the community file has to have taxon names as columns, and NO categories in the first row, so it has to be transposed. Since it cannot be transposed into a numerical matrix (it contains text as well and results in a character matrix) it has to be done stepwise: 

# transpose only the numerical parts, leaving the first column as it is
#comm2<-t(comm[,-1])
#convert it into a data frame
#comm2<-data.frame(comm2)
#assign column names to the new data frame by taking the first column from the original data frame 
#names(comm2) <- comm[,1] 


###### Serial taxa2dist
#creates a distance matrix from the aggregation file
#taxdis <- taxa2dist(agg, varstep=TRUE)

###### Parallel taxa2dist - START

  check <- TRUE
  varstep <- TRUE
  
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
#  else {
#    add <- rep(100/(ncol(agg) + check), ncol(agg) + check)
#  }
  if (!is.null(names(add))) 
    names(add) <- c("Base", names(add)[-length(add)])
  if (!check) 
    add <- c(0, add)
    
 
# Both add and x need to be distributed. But then x cannot be in the outer()   
#    addLOC <- as.matrix(add)
#    xLOC <- as.matrix(x)
    # END if (comm.rank()==0)
#} else {
#	addLOC <- NULL
#	xLOC <-NULL
#	}
    
#    add <- as.ddmatrix(addLOC)
#    x <- as.ddmatrix(xLOC)
     
    
  #out <- matrix(add[1], nrow(x), nrow(x))
  
  dim <- nrow(agg)
  
  out.d <- ddmatrix(add[1], dim, dim)
  
  
  for (i in 1:ncol(agg)) {
    #out <- out + add[i + 1] * outer(x[, i], x[, i], "!=")
    outer.d <- ddmatrix(data=add[i+1]*outer(agg[,i],agg[,i],"!="), dim, dim)
    out.d <- out.d + outer.d
  }
  

# Stop the clock
comm.print("Full program execution time for each processor",rank.source = 0)
comm.print(proc.time() - fullProgramTimer, all.rank=T)

comm.print("Size of the matrix stored in each processor (Mb)",rank.source = 0)
comm.print(object.size(out.d)/1048600,all.rank=T)

  

####### Remove this code when the data is too big to fit in processor 0
 ############# Do we want distance matrix?? 
#  taxdis <- as.matrix(out.d)
#  taxdis <- as.dist(taxdis)
#  attr(taxdis, "method") <- "taxa2dist"
#  attr(taxdis, "steps") <- add
#  attr(taxdis, "Labels") <- rownames(agg)

#comm.print("Size of the distance matrix (Mb) (notice that it is half the ddmatrix)",rank.source = 0)
#comm.print(object.size(taxdis)/1048600,rank.source = 0)


###### Parallel taxa2dist - END


#TP: The following gives the size in MBs. It is about 0.63Mb
# object.size(taxdis)/1048600

#TP: Profiling
#system.time(taxdis <- taxa2dist(agg, varstep=TRUE))
###############or
#Rprof()
#taxdis <- taxa2dist(agg, varstep=TRUE)
#Rprof(NULL)
#summaryRprof()
##############or
#Rprofmem()
#taxdis <- taxa2dist(agg, varstep=TRUE)
#Rprofmem(NULL)
#summaryRprof()
#noquote(readLines("Rprofmem.out", n = 5))


#plot(hclust(taxdis), hang = -1)



# Finish
finalize()

