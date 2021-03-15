#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(SeqArray)
sessionInfo()

argp <- arg_parser("Merge per-chromosome GDS files into single GDS file")
argp <- add_argument(argp, "--gds_file", help="chr GDS files, w/o chr identifier")
argp <- add_argument(argp, "--merged_gds_file", help="Output")
argv <- parse_args(argp)

## gds file has two parts split by chromosome identifier
gdsfile <- argv$gds_file
chr <- strsplit(c("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"), " ", fixed=TRUE)[[1]]
files <- sapply(chr, function(c) insertChromString(gdsfile, c, "gds_file"))

## write to the scratch disk of each node
gdsfile.tmp <- tempfile()
message("gds temporarily located at ", gdsfile.tmp)

## merge genotypes only (no other format or info fields)
seqMerge(files, gdsfile.tmp, fmt.var=character(), info.var=character(), storage.option="LZMA_RA")

## copy it
file.copy(gdsfile.tmp, argv$merged_gds_file)
## remove the tmp file
file.remove(gdsfile.tmp)
