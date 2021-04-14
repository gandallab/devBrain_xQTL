#! /usr/bin/env Rscript
library(argparser)
library(data.table)
library(dplyr)

p <- arg_parser("Make phenotype and sample meta input for susie")
p <- add_argument(p, "--bed", help="FastQTL bed input")
p <- add_argument(p, "--qtl_group", help="ancestry, #hcp, etc.")

args <- parse_args(p)

dat <- fread(args$bed, data.table = F)
samples <- colnames(dat)[5:ncol(dat)]
sample_meta <- data.frame("sample_id"=samples, "genotype_id"=samples, "qtl_group"=args$qtl_group)
write.table(sample_meta, paste0(dirname(args$bed),"/sample_meta.tsv"), col.names=T, row.names=F, quote=F, sep="\t")

dat2 <- dat[,c(1:4)]
dat2$phenotype_id <- dat2$ID 
dat2$group_id <- dat2$ID 
colnames(dat2)[4] <- "gene_id"
colnames(dat2)[1] <- "chromosome"
colnames(dat2)[2] <- "phenotype_pos"
dat2$strand <- rep(1, 31947)
dat3 <- dat2[,c("phenotype_id","group_id","gene_id","chromosome","phenotype_pos","strand")]

write.table(dat3, paste0(dirname(args$bed),"/phenotype_meta.tsv"), col.names=T, row.names=F, quote=F, sep="\t")
