#! /usr/bin/env Rscript
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

p <- arg_parser("Make paintor locus file")
p <- add_argument(p, "--gene", help="")
args <- parse_args(p)

# Load nominal associations
setwd("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/")
eur <- fread(paste0(args$gene, "_dir/eur_assoc.txt"), header = F, data.table = F)
amr <- fread(paste0(args$gene, "_dir/amr_assoc.txt"), header = F, data.table = F)
afr <- fread(paste0(args$gene, "_dir/afr_assoc.txt"), header = F, data.table = F)

# Intersect by common cis variants
shared <- eur %>% inner_join(amr, by = "V2") %>% inner_join(afr, by = "V2")

# Calculate zscore
# x: eur; y: amr; afr
# variants ordered by position, same as BIM
shared <- shared %>% 
    mutate(z.x = sign(V5.x)*abs(qnorm(V4.x/2))) %>% 
    mutate(z.y = sign(V5.y)*abs(qnorm(V4.y/2))) %>% 
    mutate(z = sign(V5)*abs(qnorm(V4/2)))
shared <- shared %>% select(V2, z.x, z.y, z)

# Write to locus file
colnames(shared) <- c("RSID", "ZSCORE.P1", "ZSCORE.P2", "ZSCORE.P3")
write.table(shared, paste0(args$gene, "_dir/", args$gene), col.names = T, row.names = F, sep = " ", quote = F)


id <- shared %>% select(RSID)
write.table(id, paste0(args$gene, "_dir/shared_variants.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
