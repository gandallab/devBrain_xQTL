#! /usr/bin/env Rscript

library(coloc)
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Run coloc between GWAS and ROSMAP sn eQTL")
p <- add_argument(p, "--gwas", help = "")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gwas_file", help = "")
p <- add_argument(p, "--type", help = "")

args <- parse_args(p)


test_table <- read.table(paste0(args$gwas, "_loci.tsv"), header = T)

# snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']


setwd(paste0("../out/", args$gwas, "/locus", args$locus, "/"))

# read GWAS loci 
gwas_sumstats <- fread(args$gwas_file, data.table = F)
gwas_sumstats <- gwas_sumstats %>% filter(CHROM == chr, POS > bp-500000, POS < bp+500000)

# read eqtl full assoc file
full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/ROSMAP_snRNA_eQTL/celltype-eqtl-sumstats.", args$type, ".", chr, ".tsv.gz"), data.table = F)
# Number of eQTL (significant SNP-gene pairs)
sig_assoc <- full_assoc %>% filter(significant_by_2step_FDR=='Yes')


# find genes with GWAS variants as eQTL
sig_locus <- sig_assoc %>% filter(snps %in% gwas_sumstats$ID)

# run coloc. Flip beta where neccessary
if (length(unique(sig_locus$gene_id)) > 0) {
    for(gene in unique(sig_locus$gene_id)) {
        # print(gene)
        gene_full_assoc <- full_assoc %>% filter(gene_id == gene)     
        shared <- gene_full_assoc %>% inner_join(gwas_sumstats, by = c("snps" = "ID"))
        # A1: GWAS effect; A2 GWAS other
        # REF: eQTL other; ALT: eQTL effect ???
        shared <- shared %>% filter((REF == A1 & ALT == A2) | (REF == A2 & ALT == A1))
        shared[which(shared$REF == shared$A1 & shared$ALT == shared$A2),'beta'] <- -shared[which(shared$REF == shared$A1 & shared$ALT == shared$A2),'beta']
        # create coloc data
    feature_data <- list("beta" = shared$beta,
                     "varbeta" = shared$se * shared$se,
                     "snp" = shared$snps,
                     "position" = shared$POS,
                     "type" = "quant",
                     "sdY" =1)                                   
    gwas_data <- list("beta" = shared$BETA,
                  "varbeta" = shared$SE * shared$SE,
                  "snp" = shared$snps,
                  "position" = shared$POS,
                  "type" = "quant",
                  "sdY" = 1)
    res <- coloc.abf(dataset1 = feature_data, dataset2 = gwas_data)
    if (res$summary[6] > 0.7) {
        # write.table(res$summary, paste0(gene, ".MB.coloc.sigPP4"))
        saveRDS(res, paste0(gene, ".rosmap.", args$type, ".coloc.res.rds"))

    }

    }
} else {
    # exit gracefully if no candidate target features
    quit(save = "no", status = 0)
}
