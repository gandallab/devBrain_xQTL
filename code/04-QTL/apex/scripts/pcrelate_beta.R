#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("PC-Relate")
argp <- add_argument(argp, "--out_prefix", help="")
argp <- add_argument(argp, "--gds_file", help="All chr, pruned")
argp <- add_argument(argp, "--keep_id", help="Subject IDs to keep")
argp <- add_argument(argp, "--pca_file", help="PCAiR RData")
argp <- add_argument(argp, "--n_pcs", help="Number of PC to include")
argp <- add_argument(argp, "--variant_block_size", help="Iterator size")

argv <- parse_args(argp)


# required <- c("gds_file",
#               "pca_file")
# optional <- c("n_pcs"=3,
#               "out_prefix"="pcrelate_beta",
#               "sample_include_file"=NA,
#               "variant_block_size"=1024,
#               "variant_include_file"=NA)
# config <- setConfigDefaults(config, required, optional)
# print(config)

gds <- seqOpen(argv$gds_file)
seqData <- SeqVarData(gds)

# if (!is.na(config["variant_include_file"])) {
#     filterByFile(seqData, config["variant_include_file"])
# }

pca <- getobj(argv$pca_file)
n_pcs <- min(as.integer(argv$n_pcs), length(pca$unrels))
pcs <- as.matrix(pca$vectors[,1:n_pcs])
sample.include <- samplesGdsOrder(seqData, pca$unrels)

if (!is.na(argv$keep_id)) {
    load(argv$keep_id)
    sample.id <- keep_subj_id
    sample.include <- intersect(sample.include, sample.id)
}


# create iterator
block.size <- as.integer(argv$variant_block_size)
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)

beta <- calcISAFBeta(iterator,
                     pcs=pcs,
                     sample.include=sample.include)

save(beta, file=paste0(argv$out_prefix, ".RData"))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
