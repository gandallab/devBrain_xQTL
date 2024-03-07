#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(dplyr)
library(tidyr)

p <- arg_parser("Write isoQTL BED file for FastQTL mapping with grouping")
p <- add_argument(p, "--tx2gene", help="")
p <- add_argument(p, "--bed", help="Original BED file")
p <- add_argument(p, "--info", help="GTF")
p <- add_argument(p, "--bed_out", help="BED output")

argv <- parse_args(p)

# Load data
tx2gene <- fread(argv$tx2gene, data.table = F)
bed <- fread(argv$bed, data.table = F)
info <- fread(argv$info, data.table = F)

bed <- bed %>% 
    left_join(tx2gene, by = c("ID" = "Tx")) %>% 
    unite("id", ID, Gene, sep = ":", remove = FALSE) %>% 
    left_join(info, by = c("Gene"= "V5")) %>%
    select(-2,-3,-5,-Gene,-V1)

bed <- bed[complete.cases(bed),]
# 3 genes in tx2gene but not in gencode annotation
# these tx are mapped to different genes in gencode annotation
# tx2gene generated from gentrome_decoys.tsv, which is probably from some other gencode version?
# > bed[is.na(bed$V4),1:5]
#       #Chr                                  id      1474     1496     1500
# 4989     1 ENST00000454305.1_1:ENSG00000230498  3.958807 3.803842 4.096099
# 49880   16 ENST00000565824.1_2:ENSG00000260750 -1.296910 2.946444 3.149835
# 49882   16 ENST00000563921.1_1:ENSG00000260750  1.132306 2.266702 2.177657

for (i in 1:nrow(bed)){
    if(bed[i,'V4'] == '+'){
        start <- bed[i,'V2']
        bed[i,'V2'] <- start-1
        bed[i,'V3'] <- start}
    if(bed[i,'V4'] == '-'){
        end <- bed[i,'V3']
        bed[i,'V2'] <- end-1
        bed[i,'V3'] <- end}
}
bed <- bed %>% select(1,V2,V3,id,3:(ncol(bed)-3))
colnames(bed)[1:4] <- c("#Chr", "start", "end", "ID")
bed <- bed[order(bed$'#Chr', bed$start),]
write.table(bed, argv$bed_out, col.names = T, row.names = F, quote = F, sep = "\t")