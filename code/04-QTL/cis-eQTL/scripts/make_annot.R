#! /usr/bin/env Rscript
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=T)

# Load data, subset to CHR
all <- fread(args[1], data.table = F, header = F)
coord <- fread(args[2], data.table = F, header = T)
all <- all %>% filter(V1 == args[3])
coord <- coord %>% filter(CHR == args[3])

# Write annotations to variant coordination file
# 'd' for discrete annotations
coord$CTCF_binding_site_d <- coord$enhancer_d <- coord$open_chromatin_region_d <- coord$promoter_d <- coord$promoter_flanking_region_d <- coord$TF_binding_site_d <- 0
for(i in 1:nrow(coord)) {
    temp <- all %>% filter(V4 <= coord[i,'START'], V5 >= coord[i,'START'])
    if('TF_binding_site' %in% temp$V3){
        coord[i,'TF_binding_site_d'] <- 1
    }
    if('promoter_flanking_region' %in% temp$V3){
        coord[i,'promoter_flanking_region_d'] <- 1
    }
    if('promoter' %in% temp$V3){
        coord[i,'promoter_d'] <- 1
    }
    if('open_chromatin_region' %in% temp$V3){
        coord[i,'open_chromatin_region_d'] <- 1
    }
    if('enhancer' %in% temp$V3){
        coord[i,'enhancer_d'] <- 1
    }
    if('CTCF_binding_site' %in% temp$V3){
        coord[i,'CTCF_binding_site_d'] <- 1
    }
}
coord <- coord %>% select(-c(2:4))
colnames(coord)[1] <- "SNP"
# write.table(coord, paste0(args[4], "_chr",args[3], ".tsv"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(coord, paste0(args[4], "_chr",args[3], ".tsv"), col.names = F, row.names = F, sep = "\t", quote = F)