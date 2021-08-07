#! /usr/bin/env Rscript
library(data.table)
library(dplyr)

arg <- commandArgs(trailingOnly=T)

hsq_file <- fread(paste0(arg[1], "concat_hsq.txt"), data.table = F)
# remove non-converging GCTA, NAs
hsq_file <- hsq_file %>% filter(complete.cases(hsq_file))
colnames(hsq_file) <- c("gene","h2","se","pval")
# filter for h2 and pval
sig_gene <- hsq_file %>% filter(pval < .05, h2 > 0) %>% select(gene)

sig_gene$gene <- gsub("WEIGHTS/","",sig_gene$gene)
write.table(sig_gene, paste0(arg[1], "sigificant_positive_h2.txt"), 
    col.names = F, row.names = F, quote = F, sep = "\t")

sig_gene$gene <- gsub("^", paste0(arg[1], "../WEIGHTS/"), sig_gene$gene)
sig_gene$gene <- gsub("$", ".wgt.RDat", sig_gene$gene)
write.table(sig_gene, paste0(arg[1], "../WGTLIST.txt"), 
    col.names = F, row.names = F, quote = F, sep = "\t")
