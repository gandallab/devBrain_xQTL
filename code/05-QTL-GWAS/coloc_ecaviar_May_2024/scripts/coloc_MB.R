#! /usr/bin/env Rscript

library(coloc)
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Run coloc between GWAS and metabrain eQTL")
p <- add_argument(p, "--gwas", help = "")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gwas_file", help = "")
args <- parse_args(p)


test_table <- read.table(paste0(args$gwas, "_loci.tsv"), header = T)

# snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']


setwd(paste0("../out/", args$gwas, "/locus", args$locus, "/"))

# 1. read GWAS loci 
gwas_sumstats <- fread(args$gwas_file, data.table = F)
gwas_sumstats <- gwas_sumstats %>% filter(CHROM == chr, POS > bp-500000, POS < bp+500000)

# 2. read MB full assoc file
full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-chr", chr, ".txt.gz"), data.table = F, sep = "\t")
perm <- fread("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz", data.table = F, sep = "\t")
full_assoc <- full_assoc %>% select(-c(18:20))
perm <- perm %>% select(-c(18:20))
full_assoc <- full_assoc %>% left_join(perm, by = "Gene")

# 3. filter for nominal significant eQTL
sig_assoc <- full_assoc %>% filter(MetaP.x < PvalueNominalThreshold)
sig_assoc <- sig_assoc %>% separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":")

# 4. find genes with GWAS variants as eQTL
sig_locus <- sig_assoc %>% filter(rsid %in% gwas_sumstats$ID)

# 5. run coloc. Flip beta where neccessary
if (length(unique(sig_locus$Gene)) > 0) {
    for(gene in unique(sig_locus$Gene)) {
        # print(gene)
        gene_full_assoc <- full_assoc %>% filter(Gene == gene)
        gene_full_assoc <- gene_full_assoc %>% 
        separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":") %>%
        separate(SNPAlleles.x, c("A1", "A2"), sep = "/")  
        # SNPAlleles - Alleles for the variant. Format: allele1/allele2. Note that allele1 is not always the reference allele, and allele2 is not always the effect allele. 
        rows_to_swap <- gene_full_assoc$A1 == gene_full_assoc$SNPEffectAllele.x
        temp <- gene_full_assoc$A1[rows_to_swap]
        gene_full_assoc$A1[rows_to_swap] <- gene_full_assoc$A2[rows_to_swap]
        gene_full_assoc$A2[rows_to_swap] <- temp      
        shared <- gene_full_assoc %>% inner_join(gwas_sumstats, by = c("rsid" = "ID"))
        # A1.x: MB other allele; A2.x MB effect allele
        # A1.y: GWAS effect; A2.y GWAS other
        shared <- shared %>% filter((A1.x == A2.y & A2.x == A1.y) | (A1.x == A1.y & A2.x == A2.y))
        shared[which(shared$A1.x == shared$A1.y & shared$A2.x == shared$A2.y),'MetaBeta.x'] <- -shared[which(shared$A1.x == shared$A1.y & shared$A2.x == shared$A2.y),'MetaBeta.x']
        # eqtl has rsid with A_G and G_A alleles
        dup_snp <- shared[which(duplicated(shared$rsid)), 'rsid']
        shared <- shared %>% filter(!rsid %in% dup_snp)
        # create coloc data
    feature_data <- list("beta" = shared$MetaBeta.x,
                     "varbeta" = shared$MetaSE.x * shared$MetaSE.x,
                     "snp" = shared$rsid,
                     "position" = shared$SNPPos.x,
                     "type" = "quant",
                     "sdY" =1)                                   
    gwas_data <- list("beta" = shared$BETA,
                  "varbeta" = shared$SE * shared$SE,
                  "snp" = shared$rsid,
                  "position" = shared$SNPPos.x,
                  "type" = "quant",
                  "sdY" = 1)
    res <- coloc.abf(dataset1 = feature_data, dataset2 = gwas_data)
    if (res$summary[6] > 0.7) {
        # write.table(res$summary, paste0(gene, ".MB.coloc.sigPP4"))
        saveRDS(res, paste0(gene, ".MB.coloc.res.rds"))

    }

    }
} else {
    # exit gracefully if no candidate target features
    quit(save = "no", status = 0)
}
