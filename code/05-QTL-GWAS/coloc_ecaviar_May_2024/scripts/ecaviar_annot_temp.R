library(data.table)
library(tidyverse)

res <-read.table("project-gandalm/GWAS-coloc/out/fetal_pp4.txt",header=T)
for (i in 1:nrow(res)){
test_table<-read.table(paste0("../code/",res[i,'GWAS'],"_loci.tsv"), header=T)
res[i,'CHR']<-test_table[test_table$locus==res[i,'locus'],'CHR']
res[i,'BP']<-test_table[test_table$locus==res[i,'locus'],'BP']}

res$GeneSymbol <- NA
res$GeneType <- NA

tx2gene <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/tx2gene_gencode_v33_noGeneVersion.tsv",header=T)
gtf <- fread("/u/project/gandalm/cindywen/isoform_twas/salmon/gencode.v33lift37.annotation.gene.info.tsv",data.table=F)
gtf <- gtf %>% select(ensg,V11,V12)
tx2gene <- tx2gene %>% inner_join(gtf, by = c("Gene"="ensg"))
load("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/leafviz_annot/all.introns.tested.RData")
pheno <- pheno %>% separate(ensemblID, c("ensg", "ver"), sep = "[.]")
pheno <- pheno %>% left_join(gtf, by = "ensg")
perm <- fread("/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz", data.table = F, sep = "\t")
perm <- perm %>% select(-c(18:20)) %>% select(Gene, GeneSymbol)
perm <- perm %>% separate(Gene, c("ensg", "ver"), sep = "[.]", remove = FALSE) %>% left_join(gtf, by = "ensg")


for (i in 1:nrow(res)) {
    if (res[i,'annot'] %in% c('fetal_eqtl', 'fetal_tri1_eqtl', 'fetal_tri2_eqtl')) {
        res[i,'GeneSymbol'] <- gtf[which(gtf$ensg == res[i,'gene']),'V12']
        res[i,'GeneType'] <- gtf[which(gtf$ensg == res[i,'gene']),'V11']
    } else if (res[i,'annot'] %in% c('fetal_isoqtl', 'fetal_tri1_isoqtl', 'fetal_tri2_isoqtl')) {
        res[i,'GeneSymbol'] <- tx2gene[which(tx2gene$Tx == res[i,'gene']),'V12']
        res[i,'GeneType'] <- tx2gene[which(tx2gene$Tx == res[i,'gene']),'V11']
    } else if (res[i,'annot'] %in% c('fetal_sqtl', 'fetal_tri1_sqtl', 'fetal_tri2_sqtl')) {
        res[i,'GeneSymbol'] <- pheno[which(pheno$ID == res[i,'gene']),'V12']
        res[i,'GeneType'] <- pheno[which(pheno$ID == res[i,'gene']),'V11']
    } else if (res[i,'annot'] == 'MB') {
        res[i,'GeneSymbol'] <- perm[which(perm$Gene == res[i,'gene']),'GeneSymbol']
        res[i,'GeneType'] <- perm[which(perm$Gene == res[i,'gene']),'V11']
    } 
    # else if (res[i,'annot'] == 'thistle') {
    #     res[i,'QTL'] <- 'thistle_sqtl'
    #     # probe_id <- gsub('.thistle.coloc.res.rds','',res[i,'file'])
    #     # full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/chr", res[i,'CHR'], ".txt"), data.table = F, sep = "\t")
    #     # res[i,'GeneSymbol'] <- unique(full_assoc[which(full_assoc$Probe == probe_id),'Gene'])
    #     # res[i,'GeneType'] <- gtf[which(gtf$V12 == res[i,'GeneSymbol']),'V11']
    # }
}

res_thistle <- res %>% filter(annot == "thistle")
for(chr in unique(res_thistle$CHR)) {
    full_assoc <- fread(paste0("/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/chr", chr, ".txt"), data.table = F, sep = "\t")
    for(i in which(res_thistle$CHR==chr)) {
        res_thistle[i,'GeneSymbol'] <- unique(full_assoc[which(full_assoc$Probe == res_thistle[i,'gene']),'Gene'])
        if(res_thistle[i,'GeneSymbol'] %in% gtf$V12) {
            res_thistle[i,'GeneType'] <- gtf[which(gtf$V12 == res_thistle[i,'GeneSymbol']),'V11']
        }
    }
}

res_other <- res %>% filter(annot != "thistle")

res <- rbind(res_other, res_thistle)




scz_gwas<-fread("/u/project/gandalm/shared/GWAS/SCZ.PGC3.2021/wave3_v3/PGC3_SCZ_wave3.european.autosome.public.v3.filtered.tsv.gz",data.table=F)
scz_res<-res%>%filter(GWAS=='PGC3_SCZ_wave3.european.autosome.public.v3')
scz_res<-scz_res%>%inner_join(scz_gwas, by=c("SNP_ID"="ID")) %>% select(SNP_ID,Prob_in_pCausalSet,CLPP,GWAS,locus,annot,gene,CHR,BP,GeneSymbol,GeneType,BETA,SE,PVAL)

bip_gwas<-fread("/u/project/gandalm/shared/GWAS/BIP.Mullins.2021/pgc-bip2021-all.filtered.tsv.gz",data.table=F)
bip_res<-res%>%filter(GWAS=='pgc-bip2021-all')
bip_res<-bip_res%>%inner_join(bip_gwas, by=c("SNP_ID"="ID")) %>% select(SNP_ID,Prob_in_pCausalSet,CLPP,GWAS,locus,annot,gene,CHR,BP,GeneSymbol,GeneType,BETA,SE,PVAL)

asd_gwas<-fread("/u/project/gandalm/shared/GWAS/ASD.Grove.iPSYCHPGC.2019/ASD.iPSYCHPGC.2018.filtered.tsv.gz",data.table=F)
asd_res<-res%>%filter(GWAS=='ASD.iPSYCHPGC.2018')
asd_res<-asd_res%>%inner_join(asd_gwas, by=c("SNP_ID"="ID")) %>% select(SNP_ID,Prob_in_pCausalSet,CLPP,GWAS,locus,annot,gene,CHR,BP,GeneSymbol,GeneType,BETA,SE,PVAL)

adhd_gwas<-fread("/u/project/gandalm/shared/GWAS/ADHD.Demontis.PGC.2018/ADHD.Demontis.2019.filtered.tsv.gz",data.table=F)
adhd_res<-res%>%filter(GWAS=='ADHD.Demontis.2019')
adhd_res<-adhd_res%>%inner_join(adhd_gwas, by=c("SNP_ID"="ID")) %>% select(SNP_ID,Prob_in_pCausalSet,CLPP,GWAS,locus,annot,gene,CHR,BP,GeneSymbol,GeneType,BETA,SE,PVAL)

mdd_gwas<-fread("/u/project/gandalm/shared/GWAS/MDD.Howard.PGC.2019/MDD.Howard.PGC.2019.filtered.tsv.gz",data.table=F)
mdd_res<-res%>%filter(GWAS=='MDD.Howard.PGC.2019')
mdd_res<-mdd_res%>%inner_join(mdd_gwas, by=c("SNP_ID"="ID")) %>% select(SNP_ID,Prob_in_pCausalSet,CLPP,GWAS,locus,annot,gene,CHR,BP,GeneSymbol,GeneType,BETA,SE,PVAL)

res1<-rbind(scz_res, bip_res, asd_res, adhd_res, mdd_res)
colnames(res1)[7]<-'feature'
colnames(res1)[8]<-'GWAS_loc_indexSNP_CHR'
colnames(res1)[9]<-'GWAS_loc_indexSNP_BP'
colnames(res1)[12:14]<-c("SNP_GWAS_BETA", "SNP_GWAS_SE", "SNP_GWAS_PVAL")
write.table(res1, "~/project-gandalm/GWAS-coloc/out/final_fetal_pp4_annot.tsv", col.names=T,row.names=F,quote=F,sep="\t")
