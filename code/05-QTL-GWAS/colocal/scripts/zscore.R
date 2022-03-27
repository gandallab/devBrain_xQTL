#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Extract overlapping variants, write zscore files")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gene", help = "")
p <- add_argument(p, "--trait", help = "")

args <- parse_args(p)

# 1. read input files
eqtl <- fread(paste0("/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/", args$trait, "/locus_", args$locus, "/", args$gene, "_all_pairs.txt"), data.table = F)
gwas <- read.table(paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/",
    args$trait, "/locus_",
    args$locus, "/gwas_in_bim.txt"
), header = F, stringsAsFactors = F)

# 2. overlapping variants
# > sum(unique(eqtl$V2) %in% unique(gwas$V1))
# [1] 372
# > length(unique(eqtl$V2))
# [1] 6310
eqtl_in_gwas <- eqtl %>% filter(V2 %in% gwas$V1)
# > length(unique(gwas$V1))
# [1] 556
gwas_in_eqtl <- gwas %>% filter(V1 %in% eqtl$V2)
# > dim(eqtl_in_gwas)
# [1] 372   9
# > dim(gwas_in_eqtl)
# [1] 372   9

# 3. calculate zscore
# slope/slope_se
eqtl_in_gwas <- eqtl_in_gwas %>%
    mutate(zscore = V8 / V9) %>%
    select(V2, zscore)

# this is using rs_filtered.txt sum stats, no zscore
# checking two ways to calculate zscore give the same results
# log(OR)/SE
# sign(log(OR))*abs(qnorm(P/2))
# gwas_in_eqtl <- gwas_in_eqtl %>%
# mutate(zscore1 = log(V7)/V8, zscore2 = sign(log(V7))*abs(qnorm(V9/2)))
# > cor(gwas_in_eqtl$zscore1, gwas_in_eqtl$zscore2)
# [1] 0.9999999
# now using modified sumstats.gz file, has Z, confirmed: high correlation of Z and calculated zscore
gwas_in_eqtl <- gwas_in_eqtl %>%
    # mutate(zscore = log(V7)/V8) %>%
    select(V1, V4)


# 4. sort variants
# bim <- fread("/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.12.bim", data.table = F)
# > is.unsorted(bim$V4)
# [1] FALSE
data_bim <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.bim",
    data.table = F
)
# > dim(bim)
# [1] 480110      6
# > dim(data_bim)
# [1] 8420206       6
# > is.unsorted(data_bim$V4)
# [1] TRUE
# > data_bim_12 <- data_bim %>% filter(V1 == "12")
# > is.unsorted(data_bim_12$V4)
# [1] FALSE
# make sure variants have the same coordinates
# > bim <- bim %>% inner_join(data_bim_12, by = "V2")
# > head(bim)
#   V1.x          V2      V3.x   V4.x V5.x V6.x V1.y V3.y   V4.y V5.y V6.y
# 1   12  rs77222186 0.4727766 282654    A    G   12    0 282654    A    G
# 2   12  rs76593881 0.4730730 282811    T    C   12    0 282811    T    C
# 3   12 rs114616861 0.4745155 283575    A    G   12    0 283575    A    G
# 4   12   rs2291917 0.4746370 283642    T    C   12    0 283642    T    C
# 5   12   rs7959974 0.4769184 284930    T    C   12    0 284930    T    C
# 6   12   rs7960096 0.4770823 285027    A    C   12    0 285027    A    C
# > dim(bim)
# [1] 247750     11
# > sum(bim$V4.x == bim$V4.y)
# [1] 247750

# make the variants in the same order as BIM, eQTL, and GWAS
# don't need to match with BIM actually,
# as eQTL sum stats is in coordinate order
# but gwas sum stats can be unsorted
# > is.unsorted(eqtl$V3)
# [1] FALSE
eqtl_in_gwas <- eqtl_in_gwas[order(match(eqtl_in_gwas[, 1], data_bim[, 2])), ]
gwas_in_eqtl <- gwas_in_eqtl[order(match(gwas_in_eqtl[, 1], data_bim[, 2])), ]

write.table(eqtl_in_gwas, paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/",
    args$trait, "/locus_", args$locus, "/", args$gene, "_eqtl_zscore.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

write.table(gwas_in_eqtl, paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/",
    args$trait, "/locus_", args$locus, "/", args$gene, "_gwas_zscore.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)

write.table(as.data.frame(gwas_in_eqtl$V1), paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/",
    args$trait, "/locus_", args$locus, "/", args$gene,
    "_gwas_eqtl_snp_set.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)