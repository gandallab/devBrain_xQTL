#! /usr/bin/env Rscript

# Adapted from https://github.com/francois-a/fastqtl/blob/master/R/calculateSignificanceFastQTL.R
# and https://github.com/hobrien/GENEX-FB2/blob/master/R/calulateNominalPvalueThresholds.R

# Call egene/isotx/sintron form FastQTL permutation pass
# and calculate nominal p-val threshold for qtl of each feature

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(dplyr))

p <- arg_parser("Annotates FastQTL permutation output and runs qvalue")
p <- add_argument(p, "--outdir", help="")
args <- parse_args(p)

results <- fread(paste0(args$outdir, "all.chunks.txt.gz"), data.table=F)
colnames(results) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
cat("  * FastQTL results file dimension is: ", nrow(results), " * ", ncol(results), "\n", sep="")

# 1. remove phenotypes w/o cis-variants
results <- results[complete.cases(results),]
cat("  * Complete dimension is: ", nrow(results)," * ", ncol(results), "\n", sep="")

# 2. remove duplicates
# keep phenotype with lowest nominal p-value
results <- arrange(results, npval, desc(bpval)) %>% group_by(pid) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()
cat("  * Dup-removed dimension is: ", nrow(results), " * ", ncol(results), "\n", sep="")
# Correlation between bpval and ppval
cat("  * Correlation between Beta-approximated and empirical p-values: ", round(cor(results[, 'ppval'], results[, 'bpval']), 4), "\n", sep="")

# 3. calculate q-values
# ‘signif’ rounds the values in its first argument to the specified number of significant digits
Q <- qvalue(results[,'bpval'])
results$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * Number of sig feature @ FDR 0.05: ", sum(results[, 'qval'] < 0.05), "\n", sep="")

# 4. Determine globala min(p) significance threshold and calculate nominal p-val threshold for each phenotype
ub <- sort(results[results$qval > .05, 'bpval'])[1]  # smallest p-value above FDR
lb <- -sort(-results[results$qval <= .05, 'bpval'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("  * smallest p-value above FDR: ", ub, "\n")
cat("  * largest p-value below FDR: ", lb, "\n")
cat("  * min p-value threshold @ FDR 0.05: ", pthreshold, "\n", sep="")
results[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold, results[, 'shape1'], results[, 'shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

sig_pheno <- filter(results, qval < .05)

# See OBrien code for adjusting nominal pval threshold for those genes whose best npval is larger than the threshold
# Note we only call eQTL for eGenes, so it's mediated
write.table(results, paste0(args$outdir, "all_assoc.txt"), quote=F, sep="\t", col.names=T, row.names=F)
write.table(sig_pheno, paste0(args$outdir, "sig_pheno.txt"), quote=F, sep="\t", col.names=T, row.names=F)
