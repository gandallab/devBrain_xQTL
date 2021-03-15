#! /usr/bin/env Rscript

library(argparser)
library(SeqArray)
sessionInfo()

argp <- arg_parser("Convert VCF to GDS")
argp <- add_argument(argp, "--vcf_file", help="VCF file chromosome", type="character")
argp <- add_argument(argp, "--gds_file", help="GDS file chromosome", type="character")
argp <- add_argument(argp, "--num_thread", help="Number of thread", type="character")
argv <- parse_args(argp)

seqVCF2GDS(argv$vcf_file, argv$gds_file, fmt.import="GT", storage.option="LZMA_RA", parallel=as.numeric(argv$num_thread))
