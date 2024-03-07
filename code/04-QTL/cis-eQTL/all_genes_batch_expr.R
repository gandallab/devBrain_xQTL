library(data.table)
library(tidyverse)



metadata <- read.table("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/metadata_654.tsv", header = T)
metadata <- metadata %>% filter(ancestry == "eur") %>% arrange(Age)
eqtl_samples <- read.table("~/project-gandalm/isoform_twas/eqtl_new/data/90hcp_cov_629.txt", header = T, check.names = F)
metadata <- metadata %>% filter(Subject %in% colnames(eqtl_samples))

expr <- fread("~/project-gandalm/isoform_twas/eqtl_new/data/eur/genes.280.bed.gz", data.table = F)
expr <- expr %>% select(1:4, metadata$Subject)
info <- fread("~/project-gandalm/isoform_twas/salmon/gencode.v33lift37.annotation.gene.info.tsv", data.table = F)

all_genes <- expr %>% select(ID) %>% inner_join(info, by = c("ID" = "ensg")) %>% select(ID, V12)
for(i in 1:nrow(all_genes)) {
    if (i %% 5000 == 0) {
    cat(i, "lines processed\n")
  }
    gene <- all_genes[i,'ID']
    for (j in 1:7) {
        start_id <- (j - 1)*25+1
        end_id <- start_id + 149
        if(j == 7) {
            end_id <- 280
        }
        all_genes[i, j+2] <-  rowMeans(as.matrix(expr[expr$ID == gene, (start_id+4):(end_id+4)]))[[1]]
    }
}

for(i in 1:7) {
    colnames(all_genes)[i+2] <- paste0("batch", i, "_mean_expr")
}


for(i in 1:nrow(all_genes)) {
    all_genes[i,'cor_cish2'] <- cor(as.numeric(all_genes[i,3:9]), df[,'cis_h2'])
}

colnames(all_genes)[2] <- "gene_name"
all_genes_sort <- all_genes %>% arrange(desc(abs(cor_cish2)))

info <- fread("~/project-gandalm/isoform_twas/salmon/gencode.v33lift37.annotation.gene.info.tsv", data.table = F)
info <- info %>% select(ensg, V12, V11)
all_genes_sort <- all_genes_sort %>% inner_join(info, by = c("ID" = "ensg", "gene_name" = "V12"))
colnames(all_genes_sort)[ncol(all_genes_sort)] <- "gene_type"

write.table(all_genes_sort, "~/project-gandalm/isoform_twas/eqtl_new/all_genes_batch_expr.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
