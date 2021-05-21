#! /usr/bin/env Rscript

# Refer to https://github.com/broadinstitute/gtex-pipeline/blob/e6f3c21990c199afbb9da3e4ef814274abb23ed0/qtl/leafcutter/src/map_clusters_to_genes.R

library(argparser)
library(foreach)
library(data.table)
library(dplyr)

p <- arg_parser("LeafCutter: map clusters to genes")
p <- add_argument(p, "--intron_counts_file", help="Intron counts file from LeafCutter, typically <prefix>_perind.counts.gz")
p <- add_argument(p, "--exon_file", help="File listing all unique exons in annotation. Must have columns: chr, start, end, strand, gene_id[, gene_name].")
p <- add_argument(p, "--output_name", help="Output file name")
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
argv <- parse_args(p)

exons_table<-fread(argv$exon_file, data.table = F)
intron_counts <- fread(argv$intron_counts_file, data.table = F)

introns <- intron_counts$chrom
intron_meta <- do.call(rbind, strsplit(introns,":"))
colnames(intron_meta) <- c("chr","start","end","clu")
intron_meta <- as.data.frame(intron_meta, stringsAsFactors=FALSE)
intron_meta$start <- as.numeric(intron_meta$start)
intron_meta$end <- as.numeric(intron_meta$end)

map_clusters_to_genes <- function(intron_meta, exons_table) {
  gene_df <- foreach (chr=sort(unique(intron_meta$chr)), .combine=rbind) %dopar% {

    intron_chr <- intron_meta[ intron_meta$chr==chr, ]
    exons_chr <- exons_table[exons_table$chr==chr, ]

    exons_chr$temp <- exons_chr$start
    intron_chr$temp <- intron_chr$end
    three_prime_matches <- inner_join( intron_chr, exons_chr, by="temp")

    exons_chr$temp <- exons_chr$end
    intron_chr$temp <- intron_chr$start
    five_prime_matches <- inner_join( intron_chr, exons_chr, by="temp")

    all_matches <- rbind(three_prime_matches, five_prime_matches)[ , c("clu", "gene_name")]

    all_matches <- all_matches[!duplicated(all_matches),]

    if (nrow(all_matches)==0) return(NULL)
    all_matches$clu <- paste(chr,all_matches$clu,sep=':')
    all_matches
  }

  clu_df <- gene_df %>% group_by(clu) %>% summarize(genes=paste(gene_name, collapse = ","))
  class(clu_df) <- "data.frame"
  clu_df
}

m <- map_clusters_to_genes(intron_meta, exons_table)
write.table(m, file.path(argv$output_dir, argv$output_name), sep = "\t", quote=FALSE, row.names=FALSE)
