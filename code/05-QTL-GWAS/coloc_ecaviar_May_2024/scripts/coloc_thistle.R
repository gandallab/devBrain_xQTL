#! /usr/bin/env Rscript

library(coloc)
library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Run coloc between GWAS and thistle sQTL")
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

# 2. read thistle full assoc file
full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/chr", chr, ".txt"), data.table = F, sep = "\t")

# 3. filter for nominal significant sQTL; use lenient threshold 5e-5 for thistle
sig_assoc <- full_assoc %>% filter(p < 5e-5)


# 4. find genes with GWAS variants as sQTL
sig_locus <- sig_assoc %>% filter(SNP %in% gwas_sumstats$ID)

# 5. run coloc. Flip beta where neccessary
if (length(unique(sig_locus$Probe)) > 0) {
    for(gene in unique(sig_locus$Probe)) {
        # probe can either be a gene or a leafcutter intron
        gene_full_assoc <- full_assoc %>% filter(Probe == gene)
        shared <- gene_full_assoc %>% inner_join(gwas_sumstats, by = c("SNP" = "ID"))
        shared <- shared %>% filter((A1.x == A2.y & A2.x == A1.y) | (A1.x == A1.y & A2.x == A2.y))
        # flipped alleles, flip QTL beta
        shared[which(shared$A1.x == shared$A2.y & shared$A2.x == shared$A1.y),'b'] <- -shared[which(shared$A1.x == shared$A2.y & shared$A2.x == shared$A1.y),'b']
        # creat coloc data
        feature_data <- list("beta" = shared$b,
                     "varbeta" = shared$SE.x * shared$SE.x,
                     "snp" = shared$SNP,
                     "position" = shared$BP,
                     "type" = "quant",
                     "sdY" =1)
        gwas_data <- list("beta" = shared$BETA,
                  "varbeta" = shared$SE.y * shared$SE.y,
                  "snp" = shared$SNP,
                  "position" = shared$BP,
                  "type" = "quant",
                  "sdY" = 1)
    res <- coloc.abf(dataset1 = feature_data, dataset2 = gwas_data)
    if (res$summary[6] > 0.7) {
        # write.table(res$summary, paste0(gene, ".thistle.coloc.sigPP4"))
        saveRDS(res, paste0(gene, ".thistle.coloc.res.rds"))
    }

    }
} else {
    # exit gracefully if no candidate target features
    quit(save = "no", status = 0)
}
