#! /usr/bin/env Rscript

# Edited based on Leafviz and RW's code

library(data.table)
library(argparser)
library(stringr)
library(dplyr)
library(magrittr)

p <- arg_parser("Summarise introns annotations")
p <- add_argument(p, "--annot", help="annotation output")

args <- parse_args(p)

# Load data
load(args$annot)

all.introns$verdict <- unlist(verdict.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]
all.introns$gene <- unlist(gene.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]
all.introns$ensemblID <- unlist(ensemblID.list)[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]
all.introns$transcripts <- unlist( transcripts.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]
all.introns$constitutive.score <-  unlist( constitutive.list )[ match( paste( all.introns$chr, all.introns$start, all.introns$end ), unlist(coord.list)) ]

# replace NA values/missing transcripts with "."
all.introns %<>% mutate( gene=ifelse(is.na(gene), ".", gene),
                         ensemblID=ifelse(is.na(ensemblID), ".", ensemblID),
                         transcripts=ifelse(transcripts == "", ".", transcripts),
                         constitutive.score=signif(constitutive.score, digits = 2))

# all.introns$ID <- 0
# for(i in 1:nrow(all.introns)){
#   all.introns$ID[i] <- paste0(strsplit(all.introns$chr[i],"chr")[[1]][2],":",all.introns$start[i],":",all.introns$end[i],":",all.introns$clusterID[i])
# }

save(all.introns, file=paste0(dirname(args$annot), "/all.introns.RData"))

