library(vegan)
library(foreach)
library(gRbase)
library(doMC)
library(RPostgreSQL)
library(dplyr)

x <- read.table(file="/Users/ufo/Dropbox/LW_hackaton/vegan/HMCR\ Scripts/050/aggSpecies_Percent.csv", sep = ",", quote = '"', header =T, row.names =1)



drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="test_r")
dbListConnections(drv)
dbSendQuery(con,"SET SEARCH_PATH TO vegan_test;")
dbSendQuery(con, "DROP TABLE otu_test;")

dbWriteTable(con,"otu_test",as.data.frame(x))

#Example how to create an index
dbSendQuery(con, 'CREATE INDEX ix_family ON otu_test USING btree ("Family");')
q <- fetch(dbSendQuery(con,"SELECT count(*) FROM otu_test;"))
dbDisconnect(con)

con<-src_postgres(dbname="test_r", host="localhost", user="ufo")
q <- "SELECT * FROM vegan_test.otu_test"
otu_tbl<-tbl(con, sql(q))

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



    # Create all pairs to be compared
 comb <- t(combnPrim(1:35000, 2, simplify = T))
 out <- rep(add[1], dim(comb)[1])

outerM <- function(X){
    #Define the function we will use for the outer product
    FUN <- match.fun("!=")
    #Get first part of the list
    y <- x[comb[,1],X]
    #Get second part of the list
    y1 <- x[comb[,2],X]
    #Find identical
    y2 <- FUN(y,y1)
    y2 <- add[1 + X] * y2
    return(y2)
}

s <- lapply(1:ncol(x), outerM)
dist <- Reduce('+', s)
dist <- round((out + dist), 4)
sp1 <- rownames(x)[comb[,1]]
#Get second part of the list
sp2 <- rownames(x)[comb[,2]]
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

  # Transform it into a triangular matrix, when needed:
  df <- cbind(comb,y2)
  df <- as.data.frame(df)
  zmat <- with(df, matrix(-1, ncol=max(V2), nrow=1+max(V1) ))
  zmat[with(df, cbind(V1,V2)) ] <- with(df, y2)

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

library(reshape2)

t <- cast(s4, sp1~sp2)
t <- with(s4,tapply(dist,list(sp1,sp2),"[[",1))


x <- 1:55000
m <- 2

if (length(x) == 1 && is.numeric(x))
        x <- seq(x)
    if (length(x) < m)
        stop("Error in combnPrim: n < m\n")
    NCAND <- length(x)
    NSEL <- as.integer(m)
    NSET <- as.integer(choose(NCAND, NSEL))
    ANS <- rep.int(0L, NSET * NSEL)
    res <- .C("combnC", NSEL, NCAND, NSET, ANS, DUP = TRUE, PACKAGE = "gRbase")[[4]]
    if (simplify) {
        matrix(x[res], nrow = NSEL, ncol = NSET)
    }
    else {
        res <- matrix(x[res], nrow = NSEL, ncol = NSET)
        res <- colmat2list(res)
        names(res) <- NULL
        res
    }

library(parallel)
combN <- function(X){
a <- sort(COMBS[X,])
colnames(a) <- c('v1','v2')
return(a)
}

COMBS <- mclapply(1:dim(COMBS)[1], combN, mc.cores = 4)

COMBS <- do.call("rbind", COMBS)
COMBS <- unique(COMBS)

COMBS



SPLIT <- 50000

y1 <- seq(1,SPLIT, by=1)

Y <- rep(y1, rep.int(SPLIT, SPLIT))
X <- rep(y1, times = SPLIT)

x <- seq(1,SPLIT^2, by=SPLIT+1)
y <- seq(0,SPLIT^2, by=SPLIT)
y <- y[-1]


yN <- rev(seq(SPLIT^2 - 1, SPLIT + 1, by=-(SPLIT+1)))
xN <- rev(seq((SPLIT^2)-SPLIT+1, SPLIT+1,  by=-SPLIT))

xN
yN


b <- eval(parse(text=paste("c(", paste(paste(x,':',y, sep=''),collapse = ','), ")")))
b <- eval(parse(text=paste("c(", paste(paste(xN,':',yN, sep=''),collapse = ','), ")")))

b

b1 <- seq(1,SPLIT^2)



head(x)

head(b)



if (length(b) < 10^8)
    b1 <- length(b)

if (length(b) >= 10^8)
    b1 <- seq(0, floor(length(b)/10^8) * 10^8, by=10^8)
    b1 <- c(b1, b1[length(b1)] + (length(b) - b1[length(b1)]))
    b1[1] <- 1
    b2 <- b1-1
    b2 <- b2[-1]
    b2[length(b2)] <- b2[length(b2)]+1
    b3<- cbind(b1[1:13],b2)


seri <- function(P){
X1 <- X[P[1]:P[2]]
return(X1)
}
library(parallel)
s <- mclapply(b3[1:2,],1,seri)


b1 <- b[1:10^8]



X <- X[-b]
Y <- Y[-b]

Y <- Y[b]

length(X)
head(Y, 20)
X
Y
