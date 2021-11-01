#! /usr/bin/env Rscript
library(argparser)
library(dplyr)
library(data.table)
p <- arg_parser("Make maxCPP annot file")
p <- add_argument(p, "--susie", help="")
p <- add_argument(p, "--bim", help="")
p <- add_argument(p, "--annot", help="")
args <- parse_args(p)

susie <- fread(args$susie, data.table = F)
bim <- fread(args$bim, data.table = F)
# for all variants in CS, maxCPP is the maximum CPP, or PIP here, of all eGenes
susie1 <- susie %>% 
    group_by(variant_id) %>% 
    mutate(maxCPP = max(pip)) %>% 
    select(chr, pos, variant_id, maxCPP)
susie1 <- susie1[!duplicated(susie1$variant_id),] # 122050, 4
# reference variants not in susie CS have maxCPP as 0
bim <- bim %>% left_join(susie1, by=c("V1"="chr", "V4"="pos"))
bim$maxCPP[is.na(bim$maxCPP)] <- 0
bim <- bim %>% select(V1, V4, V2, V3, maxCPP)
colnames(bim) <- c("CHR", "BP", "SNP", "CM", "maxCPP")
write.table(bim, args$annot, col.names = T, row.names = F, quote = F, sep = "\t")
