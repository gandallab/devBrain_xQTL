#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
rdata <- args[3]
# rsquared <- function(p,a){
#     r2=1-(sum((a-p)^2)/sum((a-mean(a))^2))
#     return(r2)
# }

load(rdata)
exp_genes <- gnames

# Error when no gene has cis variants in this subset (note this is not reproducible now; maybe edit later)
chr <- substring(basename(input), 4, 5)
sub <- substring(basename(input), 10, 11)
if (chr == "21" & sub == "1.") {
    out <- data.frame("")
    write.table(out, output, col.names = F, row.names = F, sep = "\t")
} else {
    dat <- read.table(input, header = T, as.is = T, check.names = F)
    goos_cor <- NULL
    for(i in 1:ncol(dat)){
        namei <- colnames(dat)[i]
        ex.flag <- which(exp_genes == namei)
        goos_cor <- c(goos_cor,cor(dat[,i],ex[,ex.flag]))
    }
    out <- data.frame(genes = colnames(dat), rsquared_goos_sva = goos_cor)
    write.table(out, output, row.names = F, quote = F, sep = "\t")
}