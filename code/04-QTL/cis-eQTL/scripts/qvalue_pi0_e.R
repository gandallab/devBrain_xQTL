#! /usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(qvalue)

############ load primary significant (discovery group)
setwd("/u/project/gandalm/cindywen/isoform_twas/")
isoqtl <- fread("isoqtl_new/results/mixed_grp_perm_70hcp/group.perm.genes.txt.gz", data.table = F)
isoqtl <- isoqtl %>% 
    separate(gene_id, c("pheno", "gene"), remove = FALSE, sep = ":") %>%
    unite("gene_qtl", gene, variant_id, remove = FALSE, sep = ":") %>%
    arrange("gene_qtl", pval_beta)
isoqtl <- isoqtl[!duplicated(isoqtl$gene_qtl),]

sqtl <- fread("sqtl_new/results/mixed_grp_perm_40hcp_1e6/group.perm.genes.txt.gz", data.table = F)
sqtl <- sqtl %>% 
    separate(gene_id, c("chr", "start", "end", "clu", "ensg", "gene_name"), sep = ":") %>%
    separate(ensg, c("gene", "version"), sep = "[.]") %>%
    unite("gene_qtl", gene, variant_id, sep = ":") %>%
    arrange("gene_qtl", pval_beta)
sqtl <- sqtl[!duplicated(sqtl$gene_qtl),]

############ load nominal association (replication group)
eqtl_nominal <- fread("eqtl_new/results/mixed_nominal_70hcp/all.chunks.txt.gz", data.table = F)
cat("Read eqtl!\n")
colnames(eqtl_nominal) <- c("pid","sid","dist","npval","slope")
eqtl_nominal <- eqtl_nominal %>% 
    unite("gene_qtl", pid, sid, remove = FALSE, sep = ":")
cat("Format nominal file done!\n")

############ qvalue
# isoQTL in eQTL
iso_e <- isoqtl %>%
    filter(qval < 0.05) %>%
    left_join(eqtl_nominal, by = "gene_qtl")
sum(is.na(iso_e$npval))
for(i in 1:nrow(iso_e)) {
    if(is.na(iso_e[i,'npval'])) {
        iso_e[i,'npval'] <- runif(1, min = 0, max = 1)
    }
}
cat("Fill in npval done!\n")

Q <- qvalue(iso_e$npval)
cat("isoqtl in eqtl pi0 is", Q$pi0, "\n")
(qvalue(iso_e$npval, lambda = seq(0.2, 0.8, 0.1)))$pi0

# sQTL in eQTL
s_e <- sqtl %>%
    filter(qval < 0.05) %>%
    left_join(eqtl_nominal, by = "gene_qtl")
sum(is.na(s_e$npval))
for(i in 1:nrow(s_e)) {
    if(is.na(s_e[i,'npval'])) {
        s_e[i,'npval'] <- runif(1, min = 0, max = 1)
    }
}
cat("Fill in npval done!\n")

Q <- qvalue(s_e$npval)
cat("sqtl in eqtl pi0 is", Q$pi0, "\n")
(qvalue(s_e$npval, lambda = seq(0.2, 0.8, 0.1)))$pi0