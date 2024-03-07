#! /usr/bin/env Rscript

# https://github.com/stephenslab/gtexresults/issues/5

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

p <- arg_parser("Calculates slope_se")
p <- add_argument(p, "--input", help="FastQTL nominal output")
p <- add_argument(p, "--output", help="")
# p <- add_argument(p, "--sig", help="Significant nominal associations")
args <- parse_args(p)

# sig <- fread(args$sig, data.table = F)
# sig <- sig %>% unite(pair, pid, sid, sep = '-', remove = FALSE)
dat <- fread(args$input, data.table = F)
colnames(dat) <- c("pid","sid","dist","npval","slope")
# dat <- dat %>% unite(pair, pid, sid, sep = '-', remove = FALSE)
# dat <- dat %>% filter(pair %in% sig$pair)
# dat$fdr <- p.adjust(dat$npval, method = 'fdr', n = as.numeric(args$num_comp))
# dat <- dat %>% filter(fdr <= 0.05)

dat$se <- NA
# GTEx fastQTL'a calculation of slope_se requires df
# tested the code below, high correlation with sebeta from GTEx
for(i in 1:nrow(dat)) {
    if(dat[i,'slope'] < 0) {
        zscore <- qnorm(dat[i,'npval']/2)
    } else {
        zscore <- qnorm(1 - dat[i,'npval']/2)
    }
    if (zscore != 0) {
        dat[i,'se'] <- dat[i,'slope']/zscore
    } else {
        dat[i,'se'] <- NA
    }
}
# torus input: loc_id, snp_id, tss_dist, pval, beata, beta_se
dat <- dat %>% select(pid, sid, dist, npval, slope, se)
write.table(dat, args$output, col.names = F, row.names = F, quote = F, sep = "\t")
