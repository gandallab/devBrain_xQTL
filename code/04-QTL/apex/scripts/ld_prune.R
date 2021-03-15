#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("LD pruning, and output GDS with subset of SNPs")
argp <- add_argument(argp, "--gds_file", help="GDS file of a chromosome")
argp <- add_argument(argp, "--keep_id", help="RData file with vector of sample.id to include")
argp <- add_argument(argp, "--genome_build", help="Genome build, used to define correlation regions")
argp <- add_argument(argp, "--maf_threshold", help="Minimum MAF for variants used in LD pruning")
argp <- add_argument(argp, "--missing_threshold", help="Maximum missing call rate for variants used in LD pruning")
argp <- add_argument(argp, "--ld_r_threshold", help="r threshold for LD pruning")
argp <- add_argument(argp, "--ld_win_size", help="Sliding window size in kb for LD pruning")
argp <- add_argument(argp, "--output_file", help="Output: GDS file with subset of pruned SNP")
argp <- add_argument(argp, "--num_thread", help="Number of thread")

argv <- parse_args(argp)

### Put these here to see what TopMed uses 
# required <- c("gds_file")
# optional <- c("autosome_only"=TRUE,
#               "exclude_pca_corr"=TRUE,
#               "genome_build"="hg38",
#               "ld_r_threshold"=0.32,
#               "ld_win_size"=10,
#               "maf_threshold"=0.01,
#               "missing_threshold"=0.01,
#               "out_file"="pruned_variants.RData",
#               "sample_include_file"=NA,
#               "variant_include_file"=NA)

# ld threshold sqrt(0.1)~=0.32


gds <- seqOpen(argv$gds_file)

# remove subjects with bad genotype quality
if (!is.na(argv$keep_id)) {
    load(argv$keep_id)
    sample.id <- keep_subj_id
    message("Using ", length(sample.id), " samples")
} else {
    sample.id <- NULL
    message("Using all samples")
}

# if (!is.na(varfile)) {
#     filterByFile(gds, varfile)
# }

# ## if we have a chromosome indicator but only one gds file, select chromosome
# if (!is.na(chr) && !bychrfile) {
#     filterByChrom(gds, chr)
# }

# filterByPass(gds)  #this is selecting 0 variants since there is no input in annotation/filter, all is NA
filterBySNV(gds) # selecting only biallelic SNV
# if (as.logical(config["exclude_pca_corr"])) {
#     filterByPCAcorr(gds, build=config["genome_build"])
# }
filterByPCAcorr(gds, build="hg19") #this appears to do nothing here? 


variant.id <- seqGetData(gds, "variant.id")
message("Using ", length(variant.id), " variants")

# auto.only <- as.logical(config["autosome_only"])
# if (chr %in% "X" & auto.only) stop("Set autosome_only=FALSE to prune X chrom variants")
maf <- as.numeric(argv$maf_threshold)
miss <- as.numeric(argv$missing_threshold)
r <- as.numeric(argv$ld_r_threshold)
win <- as.numeric(argv$ld_win_size) * 1000

set.seed(100) # make pruned SNPs reproducible
snpset <- snpgdsLDpruning(gds, sample.id=sample.id, 
                          autosome.only=TRUE, maf=maf, missing.rate=miss,
                          method="corr", slide.max.bp=win, ld.threshold=r,
                          num.thread=as.numeric(argv$num_thread))
# tried a few combos:
# function default: win=500000, r=0.2, 0.59% selected for chr21
# win=500000, r=0.32, 0.98% selected for chr21
# output are from win=50000, r=0.32 (which is sqrt(0.1)), as in Topmed pipeline, 1.52% selected 

pruned <- unlist(snpset, use.names=FALSE)
# save(pruned, file=outfile)
seqSetFilter(gds, variant.id=pruned)
seqExport(gds, argv$output_file, fmt.var=character(), info.var=character())

seqClose(gds)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
