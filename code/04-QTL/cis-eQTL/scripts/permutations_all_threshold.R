#! /usr/bin/env Rscript

# QTLtools conditional need npval threshold for all expressed genes
# refer to qtltools/script/runFDR_cis.R and our code call_perm.R
# keep genes with no cis-variants, they'll have NA as npval threshold

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(dplyr))

p <- arg_parser("npval threshold for all genes")
p <- add_argument(p, "--outdir", help="")
args <- parse_args(p)

results <- fread(paste0(args$outdir, "all.chunks.txt.gz"), data.table=F)
colnames(results) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")

MASK <- !is.na(results[,'ppval'])
Dnas <- results[!MASK,]
D <- results[MASK,]

Q <- qvalue(D[,'bpval'])
D$qval <- signif(Q$qvalues, 6)
ub <- sort(D[D$qval > .05, 'bpval'])[1]  # smallest p-value above FDR
lb <- -sort(-D[D$qval <= .05, 'bpval'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2

D[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold, D[, 'shape1'], D[, 'shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
D1 <- D[,c('pid','pval_nominal_threshold')]
D2 <- Dnas[,c('pid','ppval')]
names(D2) <- names(D1)
D3 <- rbind(D1, D2)
write.table(D3, paste0(args$outdir, "permutations_all_threshold.txt"), col.names=F, row.names=F, quote=F, sep="\t")