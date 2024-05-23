#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("Get probe in locus")
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

# read thistle full assoc file
full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/chr", chr, ".txt"), data.table = F, sep = "\t")
sig_assoc <- full_assoc %>% filter(p < 5e-5)
sig_locus <- sig_assoc %>% filter(SNP %in% gwas_sumstats$ID)


# read 1KG EUR variants, need to restrict to these to calculate LD
eur_1kg <- fread(paste0("/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.", chr, ".bim"), data.table = F)



if (length(unique(sig_locus$Probe)) > 0) {
    egene_list <- data.frame(
    'locus' = args$locus, 
    'gene' = unique(sig_locus$Probe)
    )
    
    for(gene in unique(sig_locus$Probe)) {
        gene_full_assoc <- full_assoc %>% filter(Probe == gene)
        
        shared <- gene_full_assoc %>% 
            inner_join(gwas_sumstats, by = c("SNP" = "ID")) %>% 
            inner_join(eur_1kg, by = c("SNP" = "V2"))
        
        # A1.x: QTL effect allele; A2.x: QTL other allele
        # A1.y: GWAS effect allele; A2.y: GWAS other allele
        # V5: EUR bim minor; V6: EUR bim major
        # make GWAS and QTL z respective to V5/V6
        shared <- shared %>% filter((A1.x == V5 & A2.x == V6) | (A1.x == V6 & A2.x == V5))
        shared[which(shared$A2.x == shared$V5 & shared$A1.x == shared$V6),'b'] <- -shared[which(shared$A2.x == shared$V5 & shared$A1.x == shared$V6),'b']
        shared <- shared %>% filter((A1.y == V5 & A2.y == V6) | (A1.y == V6 & A2.y == V5))
        shared[which(shared$A1.y == shared$V6 & shared$A2.y == shared$V5),'BETA'] <- -shared[which(shared$A1.y == shared$V6 & shared$A2.y == shared$V5),'BETA']
        shared <- shared %>% arrange(POS) 
        shared$qtl_z <- shared$b/shared$SE.x
        shared$gwas_z <- shared$BETA/shared$SE.y
        # eCAVIAR input files
        qtl_z <- shared %>% select(SNP, qtl_z)
        gwas_z <- shared %>% select(SNP, gwas_z)
        snps <- as.data.frame(shared$SNP)
        colnames(snps)[1] <- "V1"
        write.table(qtl_z, paste0(gene, "_thistle_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        write.table(gwas_z, paste0(gene, "_thistle_gwas_zscore.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        write.table(snps, paste0(gene, "_thistle_snps.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
    }
} else {
    egene_list <- data.frame()
}

write.table(egene_list, "locus_probe.txt", col.names = F, row.names = F, quote = F, sep = "\t")
