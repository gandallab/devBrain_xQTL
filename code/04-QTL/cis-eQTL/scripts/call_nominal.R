#! /usr/bin/env Rscript
library(argparser)
library(data.table)
library(dplyr)

p <- arg_parser("Identify FastQTL significant nominal association")
p <- add_argument(p, "--outdir", help="")
args <- parse_args(p)

results <- fread(paste0(args$outdir, "all.chunks.txt.gz"), data.table=F)
colnames(results) <- c("pid","sid","dist","npval","slope")
results$fdr <- p.adjust(results$npval, method='fdr')
write.table(results, paste0(args$outdir, "all_assoc.txt"), quote=F, sep="\t", col.names=T, row.names=F)

significant <- results %>% filter(fdr <= .05)
num_sig_feature <- length(unique(significant$pid))
df <- data.frame(num_sig_feature = c(num_sig_feature))
write.table(significant, paste0(args$outdir, "significant_assoc.txt"), quote=F, sep="\t", col.names=T, row.names=F)
write.table(df, paste0(args$outdir, "significant_feature_count.txt"), quote=F, sep="\t", col.names=F, row.names=F)
