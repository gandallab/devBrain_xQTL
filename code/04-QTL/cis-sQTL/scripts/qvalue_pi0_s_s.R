#! /usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(qvalue)

############ load primary significant (discovery group)
setwd("/u/project/gandalm/cindywen/isoform_twas/")

# sqtl_perm <- fread("sqtl_new/results/mixed_grp_perm_40hcp_1e6/group.perm.genes.txt.gz", data.table = F)
# sqtl_perm <- sqtl_perm %>%
#             filter(qval < 0.05) %>%
#             separate(group_id, c("ensg", "version"), sep = "[.]") %>%
#             unite("gene_qtl", ensg, variant_id, sep = ":") %>%
#             arrange("gene_qtl", pval_beta)

sqtl_perm <- fread("sqtl_new/results/mixed_perm_40hcp_1e6/sig_pheno_gene_info.txt", data.table = F)
sqtl_perm <- sqtl_perm %>% 
    unite("gene_qtl", ensg, sid, sep = ":") %>%
    arrange("gene_qtl", bpval)

sqtl_perm <- sqtl_perm[!duplicated(sqtl_perm$gene_qtl),]


############ load nominal association (replication group)
sqtl_nominal <- fread("sqtl_new/results/mixed_nominal_40hcp_1e6/qvalue_pi0_prep/prep.txt.gz", data.table = F, header = F)
cat("Read sqtl!\n")
colnames(sqtl_nominal) <- c("gene_qtl", "npval")

#> dim(sqtl_nominal)
#[1] 87342919        2
#length(unique(sqtl_nominal$gene_qtl))
#[1] 87342919

############ qvalue
s_s <- sqtl_perm %>% left_join(sqtl_nominal, by = "gene_qtl")
#sum(is.na(s_s$npval))
sum(is.na(s_s$npval.y))

#for(i in 1:nrow(s_s)) {
#    if(is.na(s_s[i,'npval'])) {
#        s_s[i,'npval'] <- runif(1, min = 0, max = 1)
#    }
#}
for(i in 1:nrow(s_s)) {
    if(is.na(s_s[i,'npval.y'])) {
        s_s[i,'npval.y'] <- runif(1, min = 0, max = 1)
    }
}

#Q <- qvalue(s_s$npval)
Q <- qvalue(s_s$npval.y)

cat("sqtl in sqtl pi0 is", Q$pi0, "\n")
(qvalue(s_s$npval.y, lambda = seq(0.2, 0.8, 0.1)))$pi0
