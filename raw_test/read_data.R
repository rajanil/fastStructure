library(reshape2)

d <- read.delim("testdata.ped", sep="", header=F)
sum(is.na(d))
id.cols <- c("family", "id", "pat.id", "mat.id", "sex", "phenotype")
names(d)[1:6] <- id.cols
snp.n <- (length(names(d)) - 6) / 2
if (snp.n %% 1 != 0) {
  stop("Wrong number of columns in the ped file")
}

# Name the allele columns
names(d)[-(1:6)] <- paste(c("a", "b"), rep(1:snp.n, each=2), sep=".")

d.map <- read.delim("testdata.map", sep="", header=F)
names(d.map) <- c("chromosome", "snp", "distance", "position")
if (nrow(d.map) !- snp.n) {
  stop("Wrong number of rows in the map file")
}

xi <- read.delim("testoutput_simple.3.xi", sep="", header=F)
mean.p <- read.delim("testoutput_simple.3.meanP", sep="", header=F)
mean.q <- read.delim("testoutput_simple.3.meanQ", sep="", header=F)

# Put the data set in a format that is useful to VB:
d.melt <- melt(d, id.vars=id.cols)
d.melt$location <- sub("[ab]\\.", "", d.melt$variable)
d.melt$gene <- sub("\\.[0-9]*", "", d.melt$variable)
d.shape <- dcast(d.melt, id + location ~ gene)
d.shape$g <- (d.shape$a == 2) + (d.shape$b == 2)
