#! /usr/bin/env Rscript

library(data.table)
library(argparser)

p <- arg_parser("Subset ancestries for sQTL")
p <- add_argument(p, "--combined_dir", help="cluster dir for combined ancestry data")
p <- add_argument(p, "--subj_list", help="")
p <- add_argument(p, "--outdir", help="")

args <- parse_args(p)

# Filtered, normalized, combat
mixed <- fread(paste0(args$combined_dir, "leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.tsv"), data.table=F)
rownames(mixed) <- mixed$V1
mixed <- mixed[,-1]
# No combat
noComBat <- fread(paste0(args$combined_dir, "leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj.tsv"), data.table=F)
rownames(noComBat) <- noComBat$V1
noComBat <- noComBat[,-1]
# BED
mixed.bed <- fread(paste0(args$combined_dir, "leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.bed.gz"), data.table=F)

subj.list <- read.table(args$subj_list, header = F, stringsAsFactors = F)
ancestry.in.data <- subj.list[,1] %in% colnames(mixed)
ancestry.data <- mixed[,subj.list[ancestry.in.data,]]
ancestry.noComBat <- noComBat[,subj.list[ancestry.in.data,]]
ancestry.bed <- mixed.bed[,c("#Chr","start","end","ID", subj.list[ancestry.in.data,])]

write.table(ancestry.data, paste0(args$outdir, "lc_combat.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")
write.table(ancestry.noComBat, paste0(args$outdir, "lc.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")
write.table(ancestry.bed, paste0(args$outdir, "lc_combat.bed"), col.names = T, row.names = F, quote = F, sep = "\t")
