#! /usr/bin/env Rscript
library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Filter exon GTF by mappability")
p <- add_argument(p, "--map", help="ENCODE mappability file")
p <- add_argument(p, "--exon", help="GENCODE GTF exons")
p <- add_argument(p, "--out", help="Output file")
p <- add_argument(p, "--chr", help="CHR")
argv <- parse_args(p)

# Full file
map <- fread(argv$map, data.table = F) # 295792456         4
exons <- fread(argv$exon, data.table = F) # 1379814       9

# CHR
chr <- paste0("chr", as.numeric(argv$chr))
map <- map %>% filter(V1 == chr)
exons <- exons %>% filter(V1 == chr)

# Low mappability regions: mapped to more than one region
# Only keep exons that are outside these regions
unique.map <- map %>% filter(V4 == 1)
# Merge consecutive unique map regions, if any
if(sum(unique.map$V3 %in% unique.map$V2) != 0){
    cat("Consecutive unique map regions exist!\n")
    merge <- c()
    for(i in 1:(nrow(unique.map) - 1)) {
        if (unique.map[i,3] == unique.map[i+1,2]) {
            merge <- append(merge, i)
        }
    }
    for(i in merge){
        unique.map[i+1,2] <- unique.map[i,1]
    }
    unique.map <- unique.map[-merge,]
}
# Exclude exons longer than the longest unique map region
upper.length <- range(unique.map$V3 - unique.map$V2)[2]
exons <- exons %>% mutate(length = V5-V4) %>% filter(length <= upper.length)
# Keep only exons that can be covered by unique map regions
keep <- c()
for(i in 1:nrow(exons)) {
    mapped <- unique.map %>% filter(V2 <= exons[i,'V4'], V3 >= exons[i,'V5'])
    if (nrow(mapped) >= 1) {
        keep <- append(keep, i)
    }
}
exons <- exons[keep,]
write.table(exons, argv$out, col.names = F, row.names = F, quote = F, sep = "\t")
