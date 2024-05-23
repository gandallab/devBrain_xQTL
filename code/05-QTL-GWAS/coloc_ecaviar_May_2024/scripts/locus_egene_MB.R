#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Get eGene in locus")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gwas", help = "")
p <- add_argument(p, "--gwas_file", help = "")
args <- parse_args(p)

test_table <- read.table(paste0(args$gwas, "_loci.tsv"), header = T)

# snp <- test_table[args$locus, 'SNP']
chr <- test_table[args$locus, 'CHR']
bp <- test_table[args$locus, 'BP']

setwd(paste0("../out/", args$gwas, "/locus", args$locus, "/"))

# read GWAS loci 
gwas_sumstats <- fread(args$gwas_file, data.table = F)
gwas_sumstats <- gwas_sumstats %>% filter(CHROM == chr, POS > bp-500000, POS < bp+500000)

# read MB full assoc file
full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-chr", chr, ".txt.gz"), data.table = F, sep = "\t")
perm <- fread("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz", data.table = F, sep = "\t")
full_assoc <- full_assoc %>% select(-c(18:20))
perm <- perm %>% select(-c(18:20))
full_assoc <- full_assoc %>% left_join(perm, by = "Gene")

# filter for nominal significant eQTL
sig_assoc <- full_assoc %>% filter(MetaP.x < PvalueNominalThreshold)
sig_assoc <- sig_assoc %>% separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":")

# find genes with GWAS variants as eQTL
sig_locus <- sig_assoc %>% filter(rsid %in% gwas_sumstats$ID)

# read 1KG EUR variants, need to restrict to these to calculate LD
eur_1kg <- fread(paste0("/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.", chr, ".bim"), data.table = F)



if (length(unique(sig_locus$Gene)) > 0) {
    egene_list <- data.frame(
    'locus' = args$locus, 
    'gene' = unique(sig_locus$Gene)
    )
    
    for(gene in unique(sig_locus$Gene)) {
        gene_full_assoc <- full_assoc %>% filter(Gene == gene)
        gene_full_assoc <- gene_full_assoc %>% 
        separate(SNP.x, c("snpchr", "snppos", "rsid", "snpalleles"), sep = ":") %>%
        separate(SNPAlleles.x, c("A1", "A2"), sep = "/")
        # make sure MB A1 is the ref, A2 is the effect
        rows_to_swap <- gene_full_assoc$A1 == gene_full_assoc$SNPEffectAllele.x
        temp <- gene_full_assoc$A1[rows_to_swap]
        gene_full_assoc$A1[rows_to_swap] <- gene_full_assoc$A2[rows_to_swap]
        gene_full_assoc$A2[rows_to_swap] <- temp
        # MB has multiple rows with same rsid and multiple alleles, check against GWAS and keep the entry with matching alleles
        # MB also has variants with A_G and G_A alleles, same effect allele...remove these
        shared <- gene_full_assoc %>% 
            inner_join(gwas_sumstats, by = c("rsid" = "ID")) %>% 
            inner_join(eur_1kg, by = c("rsid" = "V2"))
        dup_rsid <- shared[which(duplicated(shared$rsid)),'rsid']
        shared <- shared %>% filter(!rsid %in% dup_rsid)
        # A1.x: QTL other allele; A2.x: QTL effect allele
        # A1.y: GWAS effect allele; A2.y: GWAS other allele
        # V5: EUR bim minor; V6: EUR bim major
        # make GWAS and QTL z respective to V5/V6
        shared <- shared %>% filter((A1.x == V5 & A2.x == V6) | (A1.x == V6 & A2.x == V5))
        shared[which(shared$A1.x == shared$V5 & shared$A2.x == shared$V6),'MetaPZ.x'] <- -shared[which(shared$A1.x == shared$V5 & shared$A2.x == shared$V6),'MetaPZ.x']
        shared <- shared %>% filter((A1.y == V5 & A2.y == V6) | (A1.y == V6 & A2.y == V5))
        shared[which(shared$A1.y == shared$V6 & shared$A2.y == shared$V5),'BETA'] <- -shared[which(shared$A1.y == shared$V6 & shared$A2.y == shared$V5),'BETA']
        shared <- shared %>% arrange(snppos) # should have used GWAS hg19 POS here? match with EUR LD
        shared$Z <- shared$BETA/shared$SE
        # eCAVIAR input files
        qtl_z <- shared %>% select(rsid, MetaPZ.x)
        gwas_z <- shared %>% select(rsid, Z)
        snps <- as.data.frame(shared$rsid)
        colnames(snps)[1] <- "V1"
        write.table(qtl_z, paste0(gene, "_MB_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        write.table(gwas_z, paste0(gene, "_MB_gwas_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        write.table(snps, paste0(gene, "_MB_snps.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
    }
} else {
    egene_list <- data.frame()
}

write.table(egene_list, "locus_egene_MB.txt", col.names = F, row.names = F, quote = F, sep = "\t")
