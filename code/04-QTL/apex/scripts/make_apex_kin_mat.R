#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(dplyr)

argp <- arg_parser("Make kinship matrix for apex from GENESIS pcrelate")
argp <- add_argument(argp, "--pcrel", help="pcrelate RData")
argp <- add_argument(argp, "--out", help="Output file")
argv <- parse_args(argp)

# Note: pcrelobj$kinBtwn$kin is half of mixed_pcrelate_Matrix.RData for non-diagonal values. In Matrix.RData, diagonal values ~=1, corresponding to 0.5 in kinBtwn. Kinship matrix has 0.5 for id1 == id2. GRM has 1 for id1 == id2. So this should be specified as --kin in apex. Refer to apex documentation.
load(argv$pcrel)
mat <- pcrelobj$kinBtwn %>% select(ID1,ID2,kin)
write.table(mat, argv$out, col.names=T, row.names=F, quote=F, sep=" ")
