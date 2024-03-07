#! /usr/bin/env Rscript
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

p <- arg_parser("Merge ensembl and VEP annotations")
p <- add_argument(p, "--ensembl", help="Ensembl Regulatory Build")
p <- add_argument(p, "--vep", help="VEP output")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

# 1. Load data
ensembl <- fread(args$ensembl, data.table = F, header = F)
vep <- fread(args$vep, data.table = F, header = F)
colnames(vep)[1] <- "ID"
colnames(vep)[2] <- "Consequence"
# before I didn't use GTF and FASTA for VEP, and got 8420206 unique ID which is the total number of variants in VCF
# now added GTF and FASTA because want to specify GENCODE version, and selected only annotation with source as GENCODE GTF
# now intergenic variants do not ahve GENCODE as source, so are absent here
# > dim(vep)
# [1] 30870654        2
# > length(unique(vep$ID))
# [1] 5697720
vep <- vep %>% filter(ID %in% ensembl$V1)


# 2. Add VEP to ensembl annotations 
# follow GTEx
# note here not specifying annotation to genes, just annotate variant by any genes
# NOTE: columns names will be in reverse alphabet order
ensembl$'3_prime_UTR_variant_d' <- ensembl$'5_prime_UTR_variant_d' <- ensembl$frameshift_variant_d <- ensembl$intron_variant_d  <- ensembl$missense_variant_d <- ensembl$non_coding_transcript_exon_variant_d <- ensembl$splice_acceptor_variant_d <- ensembl$splice_donor_variant_d <- ensembl$splice_region_variant_d <- ensembl$stop_gained_d <- ensembl$synonymous_variant_d <- 0

for(i in 1:nrow(ensembl)) {
    temp <- vep %>% filter(ID == ensembl[i,1])
    vep_annot_string <- paste(temp$Consequence, collapse = ',')
    if(grepl('3_prime_UTR_variant', vep_annot_string)){
        ensembl[i,'3_prime_UTR_variant_d'] <- 1
    }
    if(grepl('5_prime_UTR_variant', vep_annot_string)){
        ensembl[i,'5_prime_UTR_variant_d'] <- 1
    }
    if(grepl('frameshift_variant', vep_annot_string)){
        ensembl[i,'frameshift_variant_d'] <- 1
    }
    if(grepl('intron_variant', vep_annot_string)){
        ensembl[i,'intron_variant_d'] <- 1
    }
    if(grepl('missense_variant', vep_annot_string)){
        ensembl[i,'missense_variant_d'] <- 1
    }
    if(grepl('non_coding_transcript_exon_variant', vep_annot_string)){
        ensembl[i,'non_coding_transcript_exon_variant_d'] <- 1
    }
    if(grepl('splice_acceptor_variant', vep_annot_string)){
        ensembl[i,'splice_acceptor_variant_d'] <- 1
    }
    if(grepl('splice_donor_variant', vep_annot_string)){
        ensembl[i,'splice_donor_variant_d'] <- 1
    }
    if(grepl('splice_region_variant', vep_annot_string)){
        ensembl[i,'splice_region_variant_d'] <- 1
    }
    if(grepl('stop_gained', vep_annot_string)){
        ensembl[i,'stop_gained_d'] <- 1
    }
    if(grepl('synonymous_variant', vep_annot_string)){
        ensembl[i,'synonymous_variant_d'] <- 1
    }
}

write.table(ensembl, args$out, col.names = F, row.names = F, sep = "\t", quote = F)
