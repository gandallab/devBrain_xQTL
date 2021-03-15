#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(GENESIS)
sessionInfo()

argp <- arg_parser("Correct PC-Relate kinship estimates, and combine sample blocks")
argp <- add_argument(argp, "--pcrelate_prefix", help="")
argp <- add_argument(argp, "--n_sample_blocks", help="")
argp <- add_argument(argp, "--sparse_threshold", help="")

argv <- parse_args(argp)

nsampblock <- as.integer(argv$n_sample_blocks)

kinSelf <- NULL
kinBtwn <- NULL
kin.thresh <- as.numeric(argv$sparse_threshold)

### correct IBD results and combine
# i know the order seems weird, but this should make sure all the correct data is loaded when needed
for (i in nsampblock:1){
    for (j in i:nsampblock){
        message('Sample Blocks ', i, ' and ', j)
        
        ## load the data
        res <- getobj(paste0(argv$pcrelate_prefix, "_block_", i, "_", j, ".RData")) 
        
        if(i == j) kinSelf <- rbind(kinSelf, res$kinSelf)
        
        # correct the IBD estimates
        res$kinBtwn <- correctK2(kinBtwn = res$kinBtwn, 
                                 kinSelf = kinSelf, 
                                 small.samp.correct = FALSE, 
                                 pcs = NULL, 
                                 sample.include = NULL)
        
        res$kinBtwn <- correctK0(kinBtwn = res$kinBtwn)
        
        # this should replace the original results, but i probably wouldn't overwrite them yet
        save(res, file=paste0(argv$pcrelate_prefix, "_block_", i, "_", j, "_corrected.RData"))

        # save results above threshold in combined file
        kinBtwn <- rbind(kinBtwn, res$kinBtwn[kin > kin.thresh])
        
        rm(res); gc()
    }
}

# save pcrelate object
pcrelobj <- list(kinSelf = kinSelf, kinBtwn = kinBtwn)
class(pcrelobj) <- "pcrelate"
save(pcrelobj, file=paste0(argv$pcrelate_prefix, ".RData"))

rm(kinBtwn, kinSelf); gc()
   
# save sparse kinship matrix
km <- pcrelateToMatrix(pcrelobj, thresh = 2*kin.thresh, scaleKin = 2)
save(km, file=paste0(argv$pcrelate_prefix, "_Matrix.RData"))

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
