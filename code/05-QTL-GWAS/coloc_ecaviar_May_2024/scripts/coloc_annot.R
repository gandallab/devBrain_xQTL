#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(argparser)

p <- arg_parser("")
p <- add_argument(p, "--gwas", help = "")
args <- parse_args(p)

res <- read.table(paste0("../out/", args$gwas, "/coloc_sigPP4.txt"))
res$GWAS <- args$gwas
res <- res %>% separate(V1, c("loc", "file"), sep="/", remove=FALSE)

# add COLOC PP4
for (i in 1:nrow(res)){
    rds <- readRDS(paste0("../out/", args$gwas, "/", res[i,'V1']))
    res[i,'PP4'] <- rds$summary[6]
}

# add GWAS locus index SNP
res$loc <- gsub('locus','',res$loc)
test_table <- read.table(paste0(args$gwas, "_loci.tsv"), header = T)
res$loc <- as.integer(res$loc)
res <- res %>% inner_join(test_table, by = c("loc"="locus"))

# add xQTL type, gene name, gene biotype
res$QTL <- NA
res$GeneSymbol <- NA
res$GeneType <- NA

tx2gene <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/tx2gene_gencode_v33_noGeneVersion.tsv",header=T)
gtf <- fread("/u/project/gandalm/cindywen/isoform_twas/salmon/gencode.v33lift37.annotation.gene.info.tsv",data.table=F)
gtf <- gtf %>% select(ensg,V11,V12)
tx2gene <- tx2gene %>% inner_join(gtf, by = c("Gene"="ensg"))
load("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/leafviz_annot/all.introns.tested.RData")
pheno <- pheno %>% separate(ensemblID, c("ensg", "ver"), sep = "[.]")
pheno <- pheno %>% left_join(gtf, by = "ensg")
perm <- fread("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz", data.table = F, sep = "\t")
perm <- perm %>% select(-c(18:20)) %>% select(Gene, GeneSymbol)
perm <- perm %>% separate(Gene, c("ensg", "ver"), sep = "[.]", remove = FALSE) %>% left_join(gtf, by = "ensg")

for (i in 1:nrow(res)) {
    if (grepl(".eqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_eqtl'
        gene_id <- gsub('.eqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- gtf[which(gtf$ensg == gene_id),'V12']
        res[i,'GeneType'] <- gtf[which(gtf$ensg == gene_id),'V11']
    } else if (grepl(".isoqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_isoqtl'
        tx_id <- gsub('.isoqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- tx2gene[which(tx2gene$Tx == tx_id),'V12']
        res[i,'GeneType'] <- tx2gene[which(tx2gene$Tx == tx_id),'V11']
    } else if (grepl(".sqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_sqtl'
        intron_id <- gsub('.sqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- pheno[which(pheno$ID == intron_id),'V12']
        res[i,'GeneType'] <- pheno[which(pheno$ID == intron_id),'V11']
    } else if (grepl("tri1_eqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_tri1_eqtl'
        gene_id <- gsub('.tri1_eqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- gtf[which(gtf$ensg == gene_id),'V12']
        res[i,'GeneType'] <- gtf[which(gtf$ensg == gene_id),'V11']
    } else if (grepl("tri2_eqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_tri2_eqtl'
        gene_id <- gsub('.tri2_eqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- gtf[which(gtf$ensg == gene_id),'V12']
        res[i,'GeneType'] <- gtf[which(gtf$ensg == gene_id),'V11']
    } else if (grepl("tri1_isoqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_tri1_isoqtl'
        tx_id <- gsub('.tri1_isoqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- tx2gene[which(tx2gene$Tx == tx_id),'V12']
        res[i,'GeneType'] <- tx2gene[which(tx2gene$Tx == tx_id),'V11']
    } else if (grepl("tri2_isoqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_tri2_isoqtl'
        tx_id <- gsub('.tri2_isoqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- tx2gene[which(tx2gene$Tx == tx_id),'V12']
        res[i,'GeneType'] <- tx2gene[which(tx2gene$Tx == tx_id),'V11']
    } else if (grepl("tri1_sqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_tri1_sqtl'
        intron_id <- gsub('.tri1_sqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- pheno[which(pheno$ID == intron_id),'V12']
        res[i,'GeneType'] <- pheno[which(pheno$ID == intron_id),'V11']
    } else if (grepl("tri2_sqtl.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'fetal_tri2_sqtl'
        intron_id <- gsub('.tri2_sqtl.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- pheno[which(pheno$ID == intron_id),'V12']
        res[i,'GeneType'] <- pheno[which(pheno$ID == intron_id),'V11']
    } else if (grepl("MB.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'MetaBrain_eqtl'
        gene_id <- gsub('.MB.coloc.res.rds','',res[i,'file'])
        res[i,'GeneSymbol'] <- perm[which(perm$Gene == gene_id),'GeneSymbol']
        res[i,'GeneType'] <- perm[which(perm$Gene == gene_id),'V11']
    } else if (grepl("thistle.coloc.res.rds", res[i,'file'], fixed = TRUE)) {
        res[i,'QTL'] <- 'thistle_sqtl'
        # probe_id <- gsub('.thistle.coloc.res.rds','',res[i,'file'])
        # full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/chr", res[i,'CHR'], ".txt"), data.table = F, sep = "\t")
        # res[i,'GeneSymbol'] <- unique(full_assoc[which(full_assoc$Probe == probe_id),'Gene'])
        # res[i,'GeneType'] <- gtf[which(gtf$V12 == res[i,'GeneSymbol']),'V11']
    }
}

res_thistle <- res %>% filter(QTL == "thistle_sqtl")
for(chr in unique(res_thistle$CHR)) {
    full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/chr", chr, ".txt"), data.table = F, sep = "\t")
    for(i in which(res_thistle$CHR==chr)) {
        probe_id <- gsub('.thistle.coloc.res.rds','',res_thistle[i,'file'])
        res_thistle[i,'GeneSymbol'] <- unique(full_assoc[which(full_assoc$Probe == probe_id),'Gene'])
        if(res_thistle[i,'GeneSymbol'] %in% gtf$V12) {
            res_thistle[i,'GeneType'] <- gtf[which(gtf$V12 == res_thistle[i,'GeneSymbol']),'V11']
        }
    }
}

res_other <- res %>% filter(QTL != "thistle_sqtl")

res <- rbind(res_other, res_thistle)

res <- res %>% select(GWAS, loc, SNP, CHR, BP, P, QTL, file, PP4, GeneSymbol, GeneType)
# res <- res %>% select(GWAS, loc, CHR, BP, QTL, file, PP4, GeneSymbol, GeneType)
write.table(res, paste0("../out/", args$gwas, "/coloc_sigPP4_annot.txt"), col.names = T, row.names = F, quote = F, sep = "\t")



