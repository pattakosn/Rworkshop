library(vegan)
library(foreach)
library(gRbase)
library(doMC)
library(RPostgreSQL)
library(dplyr)

#Read file
x <- read.table(file="/Users/ufo/Dropbox/LW_hackaton/vegan/HMCR\ Scripts/050/aggSpecies_Percent.csv", sep = ",", quote = '"', header =T, row.names =1)

#Read file
#x <- read.table(file="/megx/exchange/antonio/tmp/Rworkshop/HCMR_Scripts/050/aggSpecies_Percent.csv", sep = ",", quote = '"', header =T, row.names =1, stringsAsFactors = F)


#Create DB connection
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="test_r")

#List DB connections
dbListConnections(drv)

#Set search path to the SCHEMA
dbSendQuery(con,"SET SEARCH_PATH TO vegan_test;")

#Example of how to send a query
dbSendQuery(con, "DROP TABLE otu_test;")

#Write data to DB
dbWriteTable(con,"otu_test",as.data.frame(x))

#Example how to create an index
dbSendQuery(con, 'CREATE INDEX ix_family ON otu_test USING btree ("Family");')

#Count how many rows there are in the new created table in the DB
q <- fetch(dbSendQuery(con,"SELECT count(*) FROM otu_test;"))

#Disconnect from DB
dbDisconnect(con)

#Connect to DB using dplyr and create a tbl object
con<-src_postgres(dbname="test_r", host="localhost", user="ufo")
q <- "SELECT * FROM vegan_test.otu_test"
otu_tbl<-tbl(con, sql(q))

#Materialize tbl
x <- as.data.frame(otu_tbl, stringAsFactors=F)
rownames(x) <- x[,1]
x <- x[,-1]



varstep <- TRUE
check <- TRUE
    rich <- apply(x, 2, function(taxa) length(unique(taxa)))
    S <- nrow(x)
    if (check) {
        keep <- rich < S & rich > 1
        rich <- rich[keep]
        x <- x[, keep]
    }
    i <- rev(order(rich))
    x <- x[, i]
    rich <- rich[i]
    if (varstep) {
        add <- -diff(c(nrow(x), rich, 1))
        add <- add/c(S, rich)
        add <- add/sum(add) * 100
    }else {
        add <- rep(100/(ncol(x) + check), ncol(x) + check)
    }
    if (!is.null(names(add)))
        names(add) <- c("Base", names(add)[-length(add)])
    if (!check)
        add <- c(0, add)


combN1 <- function(Z){
    test <- function(X){
        x <- rep(1,X-1) * (Z - X ) + 1
        return(x)
    }
    test1 <- function(X){
        y <- y1[X:Z]
        return(y)
    }
    y1 <- seq(1,Z, by=1)
    s <- unlist(lapply(X=Z:1, test))
    s1 <- unlist(lapply(X=2:Z, test1))

    return(list(s,s1))
}

system.time(
comb <- combN1(90000)
)



library(gRbase)
system.time(
comb <- t(combnPrim(1:100, 2, simplify = T))
)

    # Create all pairs to be compared
# comb <- t(combnPrim(1:35000, 2, simplify = T))
# out <- rep(add[1], dim(comb)[1])
system.time(
comb <- combN1(dim(x)[1])
)
t <- combN1(5)
out <- vector("numeric",length=length(comb[[1]])) + add[1]

    b1 <- seq(0, floor(length(comb[[1]])/10^8) * 10^8, by=10^8)
    b1 <- c(b1, b1[length(b1)] + (length(comb[[1]]) - b1[length(b1)]))
    b1[1] <- 1
    b2 <- b1-1
    b2 <- b2[-1]
    b2[length(b2)] <- b2[length(b2)]+1
    b3<- cbind(b1[1:length(b1)-1],b2)

outerM <- function(X){
    outerM1 <- function(X){
#Define the function we will use for the outer product
    FUN <- match.fun("!=")
    #Get first part of the list
    y <- x[comb[[1]][b3[X,1]:b3[X,2]],1]
    #Get second part of the list
    y1 <- x[comb[[2]][b3[X,1]:b3[X,2]],1]
    #Find identical
    y2 <- FUN(y,y1)
    return(y2)
}
y2 <- unlist(lapply(1:dim(b3)[1], outerM1))

    y2 <- add[1 + X] * y2
    return(y2)
}
system.time(
s <- lapply(1:ncol(x), outerM)
    )
dist <- Reduce('+', s)
dist <- (out + dist)
sp1 <- rownames(x)[comb[[1]]]
#Get second part of the list
sp2 <- rownames(x)[comb[[2]]]
length(unique(sp1))

test <- as.data.frame(cbind(sp1,sp2,dist), stringsAsFactors = F)

library(caroline)
   #Get results
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="test_r")
dbListConnections(drv)
dbSendQuery(con,"SET SEARCH_PATH TO vegan_test;")
dbSendQuery(con, "CREATE TABLE taxa2dist_mod(sp1 text, sp2 text, dist numeric);")
dbSendQuery(con, "CREATE TABLE taxa2dist_ori(sp1 text, sp2 text, dist numeric);")


dbWriteTable(con,"taxa2dist_mod",as.data.frame(cbind(sp1, sp2, dist)), append=TRUE,overwrite=FALSE,row.names=FALSE)


library(reshape2)
s3 <- as.matrix(taxa2dist(x, varstep = T))
s4 <- melt(s3)[melt(upper.tri(s3))$value,]

s4 <- melt(s3)

colnames(s4) <- c('sp1', 'sp2', 'dist')
dbWriteTable(con,"taxa2dist_ori",s4, append=TRUE,overwrite=FALSE,row.names=FALSE)

matrix(s4,length(s4)/2,length(s4)/2)
con<-src_postgres(dbname="test_r", host="localhost", user="ufo")
q <- "SELECT * FROM vegan_test.taxa2dist"
otu_tbl<-tbl(con, sql(q))

x <- as.data.frame(otu_tbl, stringAsFactors=F)
rownames(x) <- x[,1]
x <- x[,-1]



dbDisconnect(con)


    for (i in 1:ncol(x)) {
        out <- out + add[i + 1] * outer(x[, i], x[, i], "!=")
    }
    out <- as.dist(out)
    attr(out, "method") <- "taxa2dist"
    attr(out, "steps") <- add
    if (missing(labels)) {
        attr(out, "Labels") <- rownames(x)
    } else {
        if (length(labels) != nrow(x))
            warning("Labels are wrong: needed ", nrow(x), " got ", length(labels))
        attr(out, "Labels") <- as.character(labels)
    }
    if (!check && any(out <= 0))
        warning("you used 'check=FALSE' and some distances are zero -- was this intended?")
    out

#COMBINATORICS
library(gRbase)
system.time(
comb <- t(combnPrim(1:35000, 2, simplify = T))
)


library(parallel)

combN <- function(SPLIT){
y1 <- seq(1,SPLIT, by=1)
Y <- rep(y1, rep.int(SPLIT, SPLIT))
X <- rep(y1, times = SPLIT)
x <- seq(1,SPLIT^2, by=SPLIT+1)
y <- seq(0,SPLIT^2, by=SPLIT)
y <- y[-1]


b <- eval(parse(text=paste("c(", paste(paste(x,':',y, sep=''),collapse = ','), ")")))


if (length(b) < 10^8){
    b3 <- cbind(1,length(b))
}else{
    b1 <- seq(0, floor(length(b)/10^8) * 10^8, by=10^8)
    b1 <- c(b1, b1[length(b1)] + (length(b) - b1[length(b1)]))
    b1[1] <- 1
    b2 <- b1-1
    b2 <- b2[-1]
    b2[length(b2)] <- b2[length(b2)]+1
    b3<- cbind(b1[1:length(b1)-1],b2)
}

seri <- function(P,Q){
X1 <- Q[b[b3[P,1]:b3[P,2]]]
gc()
return(X1)
}
comb <- unlist(mclapply(1:dim(b3)[1],Q=X,seri, mc.cores = 1))
comb <- cbind(unlist(mclapply(1:dim(b3)[1],Q=Y,seri, mc.cores = 10)),comb)
return(comb)
}
library(parallel)

rm(list = ls())
gc()
system.time(comb <- combN(50000)
)

SPLIT <- 5

l <- 5:1
for (i in 5:1){
cat(i)
rep(y1, rep.int(SPLIT, times = 5))
}
X <- rep(y1, times = SPLIT)




SPLIT <- 100


library(gRbase)
system.time(
comb <- t(combnPrim(1:50000, 2, simplify = T))
)


head(r)


