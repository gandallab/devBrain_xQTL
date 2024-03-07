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

eqtl <- fread("eqtl_new/results/mixed_perm_90hcp/sig_pheno.txt", data.table = F)
eqtl <- eqtl %>% 
    unite("gene_qtl", pid, sid, remove = FALSE, sep = ":")

############ load nominal association (replication group)
sqtl_nominal <- fread("sqtl_new/results/mixed_nominal_40hcp_1e6/qvalue_pi0_prep/prep.txt.gz", data.table = F, header = F)
cat("Read sqtl!\n")
colnames(sqtl_nominal) <- c("gene_qtl", "npval")

############ qvalue
# eQTL in sQTL
e_s <- eqtl %>%
    left_join(sqtl_nominal, by = "gene_qtl")
sum(is.na(e_s$npval.y))
for(i in 1:nrow(e_s)) {
    if(is.na(e_s[i,'npval.y'])) {
        e_s[i,'npval.y'] <- runif(1, min = 0, max = 1)
    }
}
cat("Fill in npval done!\n")
Q <- qvalue(e_s$npval.y)
cat("eqtl in sqtl pi0 is", Q$pi0, "\n")
(qvalue(e_s$npval.y, lambda = seq(0.2, 0.8, 0.1)))$pi0

# isoQTL in sQTL
iso_s <- isoqtl %>%
    filter(qval < 0.05) %>%
    left_join(sqtl_nominal, by = "gene_qtl")
sum(is.na(iso_s$npval))
for(i in 1:nrow(iso_s)) {
    if(is.na(iso_s[i,'npval'])) {
        iso_s[i,'npval'] <- runif(1, min = 0, max = 1)
    }
}
cat("Fill in npval done!\n")
Q <- qvalue(iso_s$npval)
cat("isoqtl in sqtl pi0 is", Q$pi0, "\n")
(qvalue(iso_s$npval, lambda = seq(0.2, 0.8, 0.1)))$pi0
