library(vegan)
library(foreach)
library(gRbase)
library(doMC)
library(RPostgreSQL)
library(dplyr)

x <- read.table(file="/Users/ufo/Dropbox/LW_hackaton/vegan/HMCR\ Scripts/010/aggSpecies_Percent.csv", sep = ",", quote = '"', header =T, row.names =1)



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



varstep <- FALSE
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
    }
    else {
        add <- rep(100/(ncol(x) + check), ncol(x) + check)
    }
    if (!is.null(names(add)))
        names(add) <- c("Base", names(add)[-length(add)])
    if (!check)
        add <- c(0, add)
    out <- matrix(add[1], nrow(x), nrow(x))

    # Create all pairs to be compared
    comb <- t(combnPrim(1:nrow(x),2, simplify = T))
    #Define the function we will use for the outer product
    FUN <- match.fun("!=")
    #Get first part of the list
    y <- x[comb[,1],1]
    #Get second part of the list
    y1 <- x[comb[,2],1]
    #Find identical
    y2 <- FUN(y,y1)
    #Get results
    y3 <- cbind(y,y1,y2)


  # Transform it into a triangular matrix, when needed:
  df <- cbind(comb,y2)
  df <- as.data.frame(df)
  zmat <- with(df, matrix(-1, ncol=max(V2), nrow=1+max(V1) ))
  zmat[with(df, cbind(V1,V2)) ] <- with(df, y2)



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

