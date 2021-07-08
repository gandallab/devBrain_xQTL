#! /usr/bin/env Rscript
library(argparser)
library(preprocessCore)

p <- arg_parser("Merge featureCounts to one rdata")
p <- add_argument(p, "--count_dir", help="Dir with featureCounts output")
p <- add_argument(p, "--gene_count", help="Number of genes in GTF")
p <- add_argument(p, "--subj_count", help="Number of subject")
p <- add_argument(p, "--out", help="Number of subject")
args <- parse_args(p)

# median.normalize.log <- function(matrix) {
#     ## normalize the sequencing data
#     log.matrix <- log2(matrix+1);

#     ##browser();
#     log.normalized.matrix <- apply(as.matrix(log.matrix), 2, function(x) { x-median(x)});

#     return(log.normalized.matrix);
# }

count_read <- function(x){
    return(length(which(x >= 1)))
}

# Read featureCounts output into matrix
file.names <- dir(args$count_dir, pattern = "count.txt$", full.names = T)

all <- matrix(0, as.numeric(args$gene_count), as.numeric(args$subj_count))
for(i in 1:length(file.names)){
    dat <- read.table(file.names[i], header = T, as.is = T , sep = "\t")
    all[,i] <- dat[,7]
}

# Remove chrX, Y, M genes
chrom <- sapply(strsplit(dat$Chr,";"),"[",1)
rm.flag <- which(chrom=="chrX" | (chrom=="chrY" | chrom=="chrM"))
if(length(rm.flag) > 0){
    all <- all[-rm.flag,]
    dat <- dat[-rm.flag,]
}

# CPM
cpm <- matrix(0, as.numeric(args$gene_count)-length(rm.flag), as.numeric(args$subj_count))
total_count <- apply(all,2,sum)
for(i in 1:ncol(all)){
    cpm[,i] <- (all[,i]/(total_count[i]/1000000))
}

# Keep genes with at least 1 CPM in at least 50% subjects
ex_count <- apply(cpm,1,count_read)
keep.flag <- which(ex_count >= as.numeric(args$subj_count)/2)

all <- cpm[keep.flag,] # all is now CPM
dat <- dat[keep.flag,]

# Convert all to TPM
gene_length <- abs(dat$Length)
# bug in original code
cpm_scale <- cpm
for (i in 1:length(gene_length)) {
    cpm_scale[i,] <- cpm_scale[i,]/gene_length[i]
}
total_count_scale <- apply(cpm_scale, 2, sum)
for(i in 1:ncol(all)){
    all[,i] <- (all[,i]/gene_length)/(total_count_scale[i]/1000000)
}

# First quantile normalize across samples, then quantile normalize across genes
all <- normalize.quantiles(as.matrix(all))
all.t <- t(all)
all2 <- normalize.quantiles(all.t)
all <- t(all2)

# Standardize
all.norm <- matrix(0, nrow(all), ncol(all))
ex_mean <- apply(all, 1, mean, na.rm=T)
for(i in 1:ncol(all)){
    all.norm[,i] <- all[,i] - ex_mean
}
ex_sd <- apply(all,1,sd,na.rm=T)
for(i in 1:ncol(all)){
    all.norm[,i] <- all.norm[,i]/ex_sd
}

# all.t: quantile normalized TPM
# ex.t: quantile normalized, standardized TPM
all.t <- t(all)
ex.t <- all.norm
ex <- t(ex.t)
gnames <- as.character(dat[,1])
save(all.t, ex, gnames, file = args$out)