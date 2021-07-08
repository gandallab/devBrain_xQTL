#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]
sub <- as.numeric(args[2])

gene_pos_file <- args[3]
genotype_dir <- args[4]
rdata <- args[5]
sva_file <- args[6]
cov_file <- args[7] #covariates, gPC, age, sex
plink_dir <- args[8]
output_file <- args[9]
sample_file <- args[10]

pos_gene <- read.table(gene_pos_file, header = T, as.is = T, sep = "\t")
allmap <- read.table(paste(genotype_dir, "chr", chr, ".bim", sep=""))


load(rdata)
exp_genes <- gnames
h2_all <- NULL
source(paste(plink_dir, "cvp.R", sep=""))

chrom_str <- paste("chr", chr, sep="")
chr_flag <- which(pos_gene$chr == chrom_str)
sub_num <- ceiling(length(chr_flag)/10)
chr_flag_min <- min(chr_flag)+(sub-1)*sub_num
chr_flag_max <- min(max(chr_flag),sub*sub_num+min(chr_flag)-1)
chr_flag <- chr_flag_min:chr_flag_max
length(chr_flag)
goos <- NULL
colgenes <- NULL
sv <- read.table(sva_file, header = T, as.is = T)
pc <- read.table(cov_file, as.is = T)
if(!is.na(pc)){
	pc <- as.matrix(pc[,])
}
if(is.na(sv[1,1])){
    sva <- pc
} else {
    sv <- as.matrix(sv)
    sva <- cbind(pc,sv)
}
for(flag_i in chr_flag){

	exp_gflag <- which(exp_genes == pos_gene$gene[flag_i])
	if(!( length(exp_gflag) == 0 )){
		cat(flag_i)
		cat("\n")
		tss1 <- max(0,pos_gene$start[flag_i]-as.numeric(args[11]))
		tss2 <- pos_gene$start[flag_i]+as.numeric(args[11])

		cis_flag <- which(allmap[,4]>tss1 & allmap[,4]<tss2)
		if(length(cis_flag)>0){
			write.table(allmap[cis_flag,c(1,4,4,4)],paste("chr",chr,"_sub",sub,".txt",sep=""),quote=F,row.names=F,col.names=F)
			system(paste(plink_dir,"plink --bfile ",genotype_dir,"chr",chr," --keep ",sample_file," --extract chr",chr,"_sub",sub,".txt --recode 12 --range --out chr",chr,"_sub",sub,sep=""))
			f_map <- paste("chr",chr,"_sub",sub,".map",sep="")
			f_ped <- paste("chr",chr,"_sub",sub,".ped",sep="")
			source(paste(plink_dir,'func_plink.R',sep=""))
			rm(ped)
			rm(map)
			source(paste(plink_dir,'func_norm.R',sep=""))
			subX=W
			M <- length(cis_flag)
			cat("Generating kinship\n")
			K <- subX %*% t(subX) / M
			cat("Saving kinship\n")
			out <- paste("chr",chr,".kinship",sep="")
			#source('func_makegrm.R')
			colgenes <- c(colgenes,pos_gene$gene[flag_i])
			phe <- ex[,exp_gflag]
			source(paste(plink_dir,'func_reml.R',sep=""))
			A <- list()
			A[[1]] <- K

			pError <- tryCatch(aiML(A=A,y=phe,verbose=F,CALC_SE=F,Var=c(0.5,0.5)),error=function(e) e)

			if(!inherits(pError,"error")){
			    reml_est <- aiML(A=A,y=phe,verbose=F,CALC_SE=F,Var=c(0.5,0.5))
			    h2g <- as.numeric(reml_est$h2)
			}else{
			    h2g=0.05
			}
			#reml_est=aiML(A=A,y=phe,verbose=F,CALC_SE=F,Var=c(0.5,0.5))

			#h2g=as.numeric(reml_est$h2)
			h2_all <- c(h2_all,h2g)
			output <- getBlups(y=phe,X=sva,Z=subX,sig2g=NA,sig2e=NA,h2=h2g)
			goos <- cbind(goos,output$gBlupLoo)
		}
	}
}
colnames(goos) <- colgenes
write.table(goos,paste(output_file),col.names=T,row.names=F,sep="\t")
out <- data.frame(gene=colgenes,h2g=h2_all)
write.table(out,paste("chr",chr,"_sub",sub,"_perm_h2g.txt",sep=""),col.names=T,row.names=F,sep="\t")