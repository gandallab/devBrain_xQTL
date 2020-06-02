#HDBR_process_metadata.R

library(data.table)
library(tidyverse)

options(stringsAsFactors = F)

meta.rnaseq = fread(("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4840/E-MTAB-4840.sdrf.txt"))

colnames(meta.rnaseq) <- colnames(meta.rnaseq) %>% str_replace("[]]", "") %>% str_replace("Characteristics", "") %>% 
  str_replace("Comment", "") %>% str_replace("Factor Value", "") %>% make.names() %>% str_replace("X.", "") %>% 
  str_replace("[.]", "_") %>% tolower()

#One row per FastQ pair --> 650 RNAseq files
meta.rnaseq <- meta.rnaseq[grepl("_1.fq.gz", meta.rnaseq$submitted_file_name),]

# 190 unique individuals with RNAseq
rnaIDs = unique(meta.rnaseq$individual)
sort(table(meta.rnaseq$organism_part))


meta.geno = fread("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4843/E-MTAB-4843.sdrf.txt")
colnames(meta.geno) <- colnames(meta.geno) %>% str_replace("[]]", "") %>% str_replace("Characteristics", "") %>% 
  str_replace("Comment", "") %>% str_replace("Factor Value", "") %>% make.names() %>% str_replace("X.", "") %>% 
  str_replace("[.]", "_") %>% tolower()

# 425 genotypes -- Keep one row per genotype run
meta.geno <- meta.geno[grepl('_Grn.idat', meta.geno$array_data.file),]

# 423 unique genotyped individuals
genoIDs = unique(meta.geno$individual)

# 183 samples with matching RNAseq and genotype data
matchedIDs = intersect(rnaIDs, genoIDs)

meta.geno <- meta.geno[meta.geno$individual %in% matchedIDs,]
meta.rnaseq <- meta.rnaseq[meta.rnaseq$individual %in% matchedIDs,]

imputedArrays = read.table('~/Downloads/imputed_HDBR.txt')[,1]
meta.geno <- meta.geno[meta.geno$assay_name %in% imputedArrays,]



obrien = c('11875','12107','12545','12546','12993','12994','13142','15240','15296','15329','15468','15533','15655','15768','16117','16286','16428','16483','16488','16548','16640','16649','16810','16826','16840','16859','16929','17013','17049','17053','17054','17068','17071','17072','17081','17087','17109','17111','17115','17130','17160','17162','17167','17175','17193','17221','17229','17264','17333','17369','17372','17475','17486','17521','17543','17629','17666','17671','17701','17753','17754','17812','17835','17921','17922','17923','17932','18015','18055','18134','18139','18153','18208','18241','18249','18266','18282','18294','18349','18355','18372','18382','18528','18529','18540','18559','18596','18611','18653','18655','18666','18687','18694','18856','18983','19031','19043','19052','1117','11396','11511','11654','11834','11844','11885','11892','11900','11907','11921','11928','11952','11971','12116','12134','12152','12153','12680','12972','13117','13144')
obrien[obrien %in% meta.geno$individual]
