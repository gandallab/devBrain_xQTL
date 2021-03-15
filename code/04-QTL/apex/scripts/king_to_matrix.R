#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Format KING results as Matrix")
argp <- add_argument(argp, "--kinship_file", help="Kinship file")
argp <- add_argument(argp, "--kinship_method", help="Kinship method")
argp <- add_argument(argp, "--out_prefix", help="Output matrix prefix")
argp <- add_argument(argp, "--sparse_threshold", help="")
argp <- add_argument(argp, "--keep_id", help="Subject ID to keep")
argv <- parse_args(argp)


# required <- c("king_file")
# optional <- c("sparse_threshold"=0.01104854, # 2^(-13/2), 5th degree
#               "out_prefix"="king_Matrix",
#               "sample_include_file"=NA,
#               "kinship_method"="king_ibdseg")
# config <- setConfigDefaults(config, required, optional)
# print(config)

# Important: because ibdseg likely does not contain all samples, as king --ibdseg is filtering sample pairs, better to have sample.id here to make output matrix contain all samples
if (!is.na(argv$keep_id)) {
    load(argv$keep_id)
    sample.id <- keep_subj_id
} else {
    sample.id <- NULL
}

if (!is.na(argv$sparse_threshold)) {
    kin.thresh <- as.numeric(argv$sparse_threshold)
} else {
    kin.thresh <- NULL
}

kin.type <- tolower(argv$kinship_method)
if (kin.type == "king_ibdseg") {
    estimator <- "PropIBD"
} else {
    estimator <- "Kinship"
}

mat <- kingToMatrix(king=argv$kinship_file,
                    estimator=estimator,
                    sample.include=sample.id,
                    thresh=kin.thresh)

if (kin.type == "king_ibdseg") {
    save(mat, file=paste0(argv$out_prefix, ".RData"))
} else {
    mat2gds(mat, paste0(argv$out_prefix, ".gds"))
}
