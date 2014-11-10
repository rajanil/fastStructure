d <- read.delim("testdata.ped", sep="", header=F)
sum(is.na(d))
names(d)[1:6] <- c("family", "id", "pat.id", "mat.id", "sex", "phenotype")
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

