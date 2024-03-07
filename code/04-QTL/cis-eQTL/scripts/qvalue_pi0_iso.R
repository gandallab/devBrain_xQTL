#! /usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(qvalue)

############ load primary significant (discovery group)
setwd("/u/project/gandalm/cindywen/isoform_twas/")
eqtl <- fread("eqtl_new/results/mixed_perm_90hcp/sig_pheno.txt", data.table = F)
eqtl <- eqtl %>% 
    unite("gene_qtl", pid, sid, remove = FALSE, sep = ":")

sqtl <- fread("sqtl_new/results/mixed_grp_perm_40hcp_1e6/group.perm.genes.txt.gz", data.table = F)
sqtl <- sqtl %>% 
    separate(gene_id, c("chr", "start", "end", "clu", "ensg", "gene_name"), sep = ":") %>%
    separate(ensg, c("gene", "version"), sep = "[.]") %>%
    unite("gene_qtl", gene, variant_id, sep = ":") %>%
    arrange("gene_qtl", pval_beta)
sqtl <- sqtl[!duplicated(sqtl$gene_qtl),]

############ load nominal association (replication group)
isoqtl_nominal <- fread("isoqtl_new/results/mixed_nominal_70hcp/all.chunks.txt.gz", data.table = F)
cat("Read isoqtl!\n")
tx2gene <- fread("salmon/tx2gene_gencode_v33_noGeneVersion.tsv", data.table = F)
colnames(isoqtl_nominal) <- c("pid","sid","dist","npval","slope")
isoqtl_nominal <- isoqtl_nominal %>% 
    left_join(tx2gene, by = c("pid" = "Tx")) %>%
    unite("gene_qtl", Gene, sid, remove = FALSE, sep = ":") %>%
    arrange("gene_qtl", npval)
isoqtl_nominal <- isoqtl_nominal[!duplicated(isoqtl_nominal$'gene_qtl'),]
cat("Format nominal file done!\n")

############ qvalue
# eQTL in isoQTL
e_iso <- eqtl %>%
    left_join(isoqtl_nominal, by = "gene_qtl")
sum(is.na(e_iso$npval.y))
for(i in 1:nrow(e_iso)) {
    if(is.na(e_iso[i,'npval.y'])) {
        e_iso[i,'npval.y'] <- runif(1, min = 0, max = 1)
    }
}
cat("Fill in npval done!\n")
Q <- qvalue(e_iso$npval.y)
cat("eqtl in isoqtl pi0 is", Q$pi0, "\n")
(qvalue(e_iso$npval.y, lambda = seq(0.2, 0.8, 0.1)))$pi0

# sQTL in isoQTL
s_iso <- sqtl %>%
    filter(qval < 0.05) %>%
    left_join(isoqtl_nominal, by = "gene_qtl")
sum(is.na(s_iso$npval))
for(i in 1:nrow(s_iso)) {
    if(is.na(s_iso[i,'npval'])) {
        s_iso[i,'npval'] <- runif(1, min = 0, max = 1)
    }
}
cat("Fill in npval done!\n")

Q <- qvalue(s_iso$npval)
cat("sqtl in isoqtl pi0 is", Q$pi0, "\n")
(qvalue(s_iso$npval, lambda = seq(0.2, 0.8, 0.1)))$pi0
