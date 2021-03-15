#! /usr/bin/env Rscript
library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(GENESIS)
sessionInfo()

argp <- arg_parser("PC-Relate")
argp <- add_argument(argp, "--gds_file", help="All chr, pruned")
argp <- add_argument(argp, "--beta_file", help="ISAF beta")
argp <- add_argument(argp, "--pca_file", help="pcair RData")
argp <- add_argument(argp, "--keep_id", help="")
argp <- add_argument(argp, "--n_pcs", help="")
argp <- add_argument(argp, "--out_prefix", help="")
argp <- add_argument(argp, "--n_sample_blocks", help="Number of sample blocks")
argp <- add_argument(argp, "--variant_block_size", help="")
# argp <- add_argument(argp, "--segment",help="Calculated from sample blocks (???)")
argp <- add_argument(argp, "--block1", help="")
argp <- add_argument(argp, "--block2", help="")
argv <- parse_args(argp)
# cat(">>> TopmedPipeline version ", argv$version, "\n")
# config <- readConfig(argv$config)
# segment <- argv$segment

# required <- c("gds_file",
#               "pca_file",
#               "beta_file")
# optional <- c("n_pcs"=3,
#               "out_prefix"="pcrelate",
#               "n_sample_blocks"=1,
#               #"sample_block_size"=10000,
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
sample.include <- samplesGdsOrder(seqData, rownames(pcs))

if (!is.na(argv$keep_id)) {
    load(argv$keep_id)
    sample.id <- keep_subj_id
    sample.include <- intersect(sample.include, sample.id)
}

# load betas
beta <- getobj(argv$beta_file)

# create iterator
block.size <- as.integer(argv$variant_block_size)
iterator <- SeqVarBlockIterator(seqData, variantBlock=block.size)

# create sample blocks
#nsampblock <- ceiling(length(sample.include)/sample.block.size)
nsampblock <- as.integer(argv$n_sample_blocks)
# if (nsampblock > 1) {
#     samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), nsampblock)))

#     # map segment number to sample block numbers
#     ## jobs <- list()
#     ## s <- 1
#     ## for (i in 1:nsampblock) {
#     ##     for (j in i:nsampblock) {
#     ##         jobs[[s]] <- c(i,j)
#     ##         s <- s + 1
#     ##     }
#     ## }
#     jobs <- c(combn(1:nsampblock, 2, simplify=FALSE),
#               lapply(1:nsampblock, function(x) c(x,x)))
#     i <- jobs[[segment]][1]
#     j <- jobs[[segment]][2]
# } else {
#     samp.blocks <- list(sample.include)
#     i <- j <- 1
# }

if (nsampblock > 1) {
    samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), nsampblock)))
    i <- as.numeric(argv$block1)
    j <- as.numeric(argv$block2)
} else {
    samp.blocks <- list(sample.include)
    i <- j <- 1
}

out <- pcrelateSampBlock(iterator,
                         betaobj=beta,
                         pcs=pcs,
                         sample.include.block1=samp.blocks[[i]],
                         sample.include.block2=samp.blocks[[j]])

save(out, file=paste0(argv$out_prefix, "_block_", i, "_", j, ".RData"))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
