#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
count_dir <- args[2]
out_dir <- args[3]

load(input)
g <- data.frame(gene = gnames, order = 1:length(gnames))
allfiles <- list.files(count_dir, pattern = ".count.txt", full.names = T)
dat <- read.table(allfiles[1], header = T, sep = "\t")
all <- merge(g, dat, by.x = 1, by.y = 1)
all <- all[order(all$order),]
chrom <- sapply(strsplit(as.character(all$Chr),";"),"[",1)
start <- sapply(strsplit(as.character(all$Start),";"),"[",1)
end <- sapply(strsplit(as.character(all$End),";"),tail,1)

out <- data.frame(gene = all$gene, chr = chrom, start = start, end = end)
write.table(out, paste(out_dir, "gene_pos.txt", sep=""), row.names = F, sep = "\t", quote = F)