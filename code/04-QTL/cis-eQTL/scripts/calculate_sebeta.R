#! /usr/bin/env Rscript

# https://github.com/stephenslab/gtexresults/issues/5

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

p <- arg_parser("Calculates slope_se")
p <- add_argument(p, "--input", help="FastQTL nominal output")
p <- add_argument(p, "--output", help="")
args <- parse_args(p)

dat <- fread(args$input, data.table = F)
dat$se <- NA
# GTEx fastQTL'a calculation of slope_se requires df
# tested the following code, with and without slope>0, <0 difference, in GTEx v8 Cortex egene file
# differentiating slope>0, <0 gives high correlation with GTEx slope_se (>0.99)
# no differentiating cor~0.03
for(i in 1:nrow(dat)) {
    if(dat[i,5] < 0) {
        zscore <- qnorm(dat[i,4]/2)
    } else {
        zscore <- qnorm(1 - dat[i,4]/2)
    }
    if (zscore != 0) {
        dat[i,'se'] <- dat[i,5]/zscore
    } else {
        dat[i,'se'] <- NA
    }
}
write.table(dat, args$output, col.names = F, row.names = F, quote = F, sep = "\t")
