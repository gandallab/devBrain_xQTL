#! /usr/bin/env Rscript
library(data.table)
library(dplyr)

arg <- commandArgs(trailingOnly = T)

weights <- fread(arg[1], data.table = F, header = F)
bed <- fread(arg[2], data.table = F, header = T)[, 2:4]
sig_gene <- fread(arg[3], data.table = F, header = F)
gtf_gene <- fread(arg[4], data.table = F)[, c(10, 13)]

bed <- bed %>% filter(Gene_Symbol %in% sig_gene$V1)
gtf_gene <- gtf_gene %>% filter(ensg %in% sig_gene$V1)

weights$gene <- NA
for (i in 1:nrow(weights)) {
    weights[i, "gene"] <- strsplit(basename(weights[i, 1]), split = "[.]")[[1]][1]
}
df <- weights %>%
    left_join(bed, c("gene" = "Gene_Symbol")) %>%
    left_join(gtf_gene, c("gene" = "ensg"))

colnames(df) <- c("WGT", "gene", "CHR", "P0", "ID")
df$P1 <- df$P0 + 1
df$PANEL <- arg[6]
df <- df %>% select(PANEL, WGT, ID, CHR, P0, P1)
df$WGT <- basename(df$WGT)

write.table(df, arg[5], col.names = T, row.names = F, quote = F, sep = "\t")