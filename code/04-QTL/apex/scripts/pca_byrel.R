#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("PCA on unrelated set, project relatives")
argp <- add_argument(argp, "--rel_file", help="partitioned related subset")
argp <- add_argument(argp, "--unrel_file", help="partitioned unrelated subset")
argp <- add_argument(argp, "--gds_file", help="All chr GDS file")
argp <- add_argument(argp, "--keep_id", help="Subject ID to keep")
argp <- add_argument(argp, "--out_file", help="Output PCA")
argp <- add_argument(argp, "--out_file_unrel", help="Output PCA unrelated")
argp <- add_argument(argp, "--n_pcs", help="Number of PCs to return")
argp <- add_argument(argp, "--num_thread", help="")

argv <- parse_args(argp)

# required <- c("gds_file",
#               "related_file",
#               "unrelated_file")
# optional <- c("n_pcs"=32,
#               "out_file"="pca.RData",
#               "out_file_unrel"="pca_unrel.RData",
#               "sample_include_file"=NA,
#               "variant_include_file"=NA)
# config <- setConfigDefaults(config, required, optional)
# print(config)

gds <- seqOpen(argv$gds_file)

rels <- getobj(argv$rel_file)
unrels <- getobj(argv$unrel_file)

if (!is.na(argv$keep_id)) {
    load(argv$keep_id)
    sample.id <- keep_subj_id
    rels <- intersect(rels, sample.id)
    unrels <- intersect(unrels, sample.id)
}
message("Using ", length(unrels), " unrelated and ", length(rels), " related samples")

# if (!is.na(config["variant_include_file"])) {
#     variant.id <- getobj(config["variant_include_file"])
#     message("Using ", length(variant.id), " variants")
# } else {
#     variant.id <- NULL
# }

# number of PCs to return
n_pcs <- min(as.integer(argv$n_pcs), length(unrels))

# run PCA on unrelated set
message("PCA on unrelated set")
# nt <- countThreads()
pca.unrel <- snpgdsPCA(gds, sample.id=unrels, algorithm="randomized", eigen.cnt=n_pcs, num.thread=as.numeric(argv$num_thread))
save(pca.unrel, file=argv$out_file_unrel)

# project values for relatives
message("PCA projection for related set")
snp.load <- snpgdsPCASNPLoading(pca.unrel, gdsobj=gds, num.thread=as.numeric(argv$num_thread))
samp.load <- snpgdsPCASampLoading(snp.load, gdsobj=gds, sample.id=rels, num.thread=as.numeric(argv$num_thread))

# combine unrelated and related PCs and order as in GDS file
eigenvect <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
rownames(eigenvect) <- c(pca.unrel$sample.id, samp.load$sample.id)
seqSetFilter(gds, sample.id=rownames(eigenvect), verbose=FALSE)
sample.id <- seqGetData(gds, "sample.id")
samp.ord <- match(sample.id, rownames(eigenvect))
eigenvect <- eigenvect[samp.ord,]

# output object
pca <- list(vectors=eigenvect,
            values=pca.unrel$eigenval[1:n_pcs],
            varprop=pca.unrel$varprop[1:n_pcs],
            rels=rels,
            unrels=unrels)

save(pca, file=argv$out_file)

seqClose(gds)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
