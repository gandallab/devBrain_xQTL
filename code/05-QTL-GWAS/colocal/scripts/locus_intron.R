#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Get intron in locus")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--trait", help = "")
args <- parse_args(p)

locus_file <- read.table(paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/",
    args$trait, "/locus_", args$locus, "/gwas_in_bim.txt"
),
header = F, stringsAsFactors = F
)
# > head(locus_File)
#           V1     V2 V3      V4 V5 V6      V7        V8        V9
# 1  rs1006690 0.9612 12 2375636  C  T 0.97054 0.0223380 1.808e-01
# 2  rs1006737 0.6839 12 2345295  G  A 0.91583 0.0101120 3.470e-18
# 3 rs10160993 0.5467 12 2499565  C  G 0.96803 0.0094896 6.160e-04
# 4  rs1016388 0.5924 12 2321868  A  T 0.92744 0.0095085 2.340e-15
# 5  rs1024582 0.3290 12 2402246  A  G 1.09210 0.0101640 4.490e-18
# 6 rs10466899 0.9920 12 2482943  C  T 0.95823 0.0357000 2.320e-01
egene <- fread(
    "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/significant_assoc.txt",
    data.table = F
)
egene_in_locus <- egene %>% filter(sid %in% locus_file$V1)
egene_list <- as.data.frame(unique(egene_in_locus$pid))
write.table(egene_list, paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_sqtl/",
    args$trait, "/locus_", args$locus, "/locus_intron.txt"
),
col.names = F, row.names = F, quote = F, sep = "\t"
)