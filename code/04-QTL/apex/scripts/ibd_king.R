#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("IBD with KING-robust")
argp <- add_argument(argp, "--gds_file", help="LD pruned GDS file, all chr")
argp <- add_argument(argp, "--out_file", help="Output KING-robust GDS file")
argp <- add_argument(argp, "--keep_id", help="Subject ID to keep")
argp <- add_argument(argp, "--num_thread", help="Number of threads to use")
argv <- parse_args(argp)


# required <- c("gds_file")
# optional <- c("out_file"="ibd_king.gds",
#               #"out_file"="ibd_king.RData",
#               "sample_include_file"=NA,
#               "variant_include_file"=NA)
# config <- setConfigDefaults(config, required, optional)
# print(config)

gds <- seqOpen(argv$gds_file)

if (!is.na(argv$keep_id)) {
    load(argv$keep_id)
    sample.id <- keep_subj_id
} else {
    sample.id <- NULL
    message("Using all samples")
}

# if (!is.na(config["variant_include_file"])) {
#     variant.id <- getobj(config["variant_include_file"])
#     message("Using ", length(variant.id), " variants")
# } else {
#     variant.id <- NULL
#     message("Using all variants")
# }

ibd <- snpgdsIBDKING(gds, sample.id=sample.id, snp.id=NULL,
                     num.thread=as.numeric(argv$num_thread))

#save(ibd, file=config["out_file"])
list2gds(ibd, argv$out_file)

seqClose(gds)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
