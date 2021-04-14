#! /usr/bin/env Rscript
library(data.table)
library(dplyr)

arg <- commandArgs(trailingOnly=T)

weights <- fread(arg[1], data.table=F, hea=F)
bed <- fread(arg[2], data.table=F, hea=T)[,2:4]
sig_gene <- fread(arg[3], data.table=F, hea=F)
gtf_gene <- fread(arg[4], data.table=F)[,c(10,13)]

bed <- bed %>% filter(Gene_Symbol %in% sig_gene$V1)
gtf_gene <- gtf_gene %>% filter(ensg %in% sig_gene$V1)

weights$gene <- substring(weights$V1, 64, 78)
df <- weights %>% left_join(bed, c("gene"="Gene_Symbol"))
df2 <- df %>% left_join(gtf_gene, c("gene"="ensg"))

colnames(df2) <- c("WGT","gene","CHR","P0","ID")
df2$P1 <- df2$P0 + 1
df2 <- df2[,c("WGT","ID","CHR","P0","P1")]
df2$WGT <- basename(df2$WGT)

write.table(df2, arg[5], col.names=T, row.names=F, quote=F, sep="\t")
