#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Get feature in locus")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gwas", help = "")
p <- add_argument(p, "--gwas_file", help = "")
p <- add_argument(p, "--annot", help = "fetal xQTL")

args <- parse_args(p)

test_table <- read.table(paste0(args$gwas, "_loci.tsv"), header = T)

# snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']

setwd(paste0("../out_1MB/", args$gwas, "/locus", args$locus, "/"))

# read GWAS loci 
gwas_sumstats <- fread(args$gwas_file, data.table = F)
gwas_sumstats <- gwas_sumstats %>% filter(CHROM == chr, POS > bp-1000000, POS < bp+1000000)

# read xQTL
if (args$annot == "eqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "isoqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "sqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/gtex_chunk", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri1_eqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/tri1_25_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T1-nominal_significant_assoc.txt", data.table = F)
} else if (args$annot == "tri2_eqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/tri2_15_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T2-nominal_significant_assoc.txt", data.table = F)
} else if (args$annot == "tri1_isoqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri1_nominal_35hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri1_nominal_35hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri2_isoqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri2_nominal_20hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri2_nominal_20hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri1_sqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri1_nominal_15hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri1_nominal_15hcp/significant_assoc.txt", data.table = F)
} else if (args$annot == "tri2_sqtl") {
    full_assoc <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri2_nominal_10hcp/gtex_chr", chr, ".txt.gz"), data.table = F, sep = "\t")
    full_sig <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri2_nominal_10hcp/significant_assoc.txt", data.table = F)
}

colnames(full_assoc) <- c("gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se")
# trimester files has NA
full_assoc <- full_assoc[complete.cases(full_assoc),]

# get A1/A2 of QTL sum stats
if (args$annot %in% c("eqtl", "isoqtl", "sqtl")) {
    qtl_bim <- fread("/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test.bim", data.table = F)
} else if (grepl("tri", args$annot)) {
    qtl_bim <- fread("/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.bim", data.table = F)
}
full_assoc <- full_assoc %>% inner_join(qtl_bim, by = c("variant_id" = "V2"))

full_assoc <- full_assoc %>% unite("pair", gene_id, variant_id, sep = "-", remove = FALSE)
full_sig <- full_sig %>% unite("pair", pid, sid, sep = "-", remove = FALSE)
sig_assoc <- full_assoc %>% filter(pair %in% full_sig$pair)



sig_locus <- sig_assoc %>% filter(variant_id %in% gwas_sumstats$ID)


# read 1KG EUR variants, need to restrict to these to calculate LD for GWAS
eur_1kg <- fread(paste0("/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.", chr, ".bim"), data.table = F)



if (length(unique(sig_locus$gene_id)) > 0) {
    egene_list <- data.frame(
    'locus' = args$locus, 
    'gene' = unique(sig_locus$gene_id)
    )
    
    for(gene_i in unique(sig_locus$gene_id)) {
        gene_full_assoc <- full_assoc %>% filter(gene_id == gene_i)
        
        shared <- gene_full_assoc %>% 
            inner_join(gwas_sumstats, by = c("variant_id" = "ID")) %>% 
            inner_join(eur_1kg, by = c("variant_id" = "V2"))
        
        # make sure there are shared variants after joining with 1kg
        if (nrow(shared > 0)) {
            # V5.x: QTL effect allele; V6.x: QTL other allele
            # A1: GWAS effect allele; A2: GWAS other allele
            # V5.y: EUR bim minor; V6.y: EUR bim major
            # fetal trimester: make GWAS and QTL z respective to V5.y/V6.y
            # fetal e/iso/sQTL: QTL z respective to V5.x/V6.x
            shared <- shared %>% filter((V5.x == V5.y & V6.x == V6.y) | (V5.x == V6.y & V6.x == V5.y))
            if (grepl("tri", args$annot)) {
                shared[which(shared$V5.x == shared$V6.y & shared$V6.x == shared$V5.y),'slope'] <- -shared[which(shared$V5.x == shared$V6.y & shared$V6.x == shared$V5.y),'slope']
            } 
            shared <- shared %>% filter((A1 == V5.y & A2 == V6.y) | (A1 == V6.y & A2 == V5.y))
            shared[which(shared$A1 == shared$V6.y & shared$A2 == shared$V5.y),'BETA'] <- -shared[which(shared$A1 == shared$V6.y & shared$A2 == shared$V5.y),'BETA']
            shared <- shared %>% arrange(POS) 
            shared$qtl_z <- shared$slope/shared$slope_se
            shared$gwas_z <- shared$BETA/shared$SE
            # eCAVIAR input files
            qtl_z <- shared %>% select(variant_id, qtl_z)
            gwas_z <- shared %>% select(variant_id, gwas_z)
            snps <- as.data.frame(shared$variant_id)
            colnames(snps)[1] <- "V1"
            write.table(qtl_z, paste0(gene_i, "_fetal_", args$annot, "_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
            write.table(gwas_z, paste0(gene_i, "_fetal_", args$annot, "_gwas_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
            write.table(snps, paste0(gene_i, "_fetal_", args$annot, "_snps.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        } else {
            egene_list <- egene_list %>% filter(gene != gene_i)
        }
        
    }
} else {
    egene_list <- data.frame()
}

write.table(egene_list, paste0("locus_fetal_", args$annot, ".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
