#! /usr/bin/env Rscript
# library(phenix)
library(PHENIX)

args <- commandArgs(trailingOnly=T)
chr <- args[1]
sub <- args[2]
pr_file <- args[3]
output <- args[4]
rdata <- args[5]
cvblup_file <- args[6]
tool_dir <- args[7]
trans_dir <- args[8]
load(rdata)
source(paste(tool_dir,'make_smartsva.R',sep=""))
exp_genes <- gnames

g <- read.table(paste(cvblup_file, sep = ""), header = T, check.names = F)
pr <- read.table(pr_file, as.is = T, header = T, sep = "\t")
h2 <- read.table(paste(trans_dir, "chr",chr,"_sub",sub,"_perm_h2g.txt",sep = ""), as.is = T, header = T, sep = "\t")
coln <- NULL
allcor <- NULL
allcoef <- NULL
for(i in 1:ncol(g)){
	if(!is.na(g[1,i]) & (pr[i,2]>0.1 &(h2[i,2]>0 & h2[i,2]<1)){
		tryCatch({
		g_sva=make_sva(ex,100,g[,i])
		g_norm=(g[,i]-mean(g[,i]))/sd(g[,i])
		if(!is.na(g_sva[1,1])){
	        rlm=lm(ex~g_sva[,1:20])
	        rQN=quantnorm(as.matrix(resid(rlm)))} else{rQN=ex}
	        test=lm(rQN~g_norm)
	        allcor=cbind(allcor,sapply(summary(test),function(x) x$coefficients[2,4]))
	        allcoef=cbind(allcoef,sapply(summary(test),function(x) x$coefficients[2,1]))
	        coln=c(coln,colnames(g)[i])
		}, error=function(e){})
		#test=apply(ex_norm,2,cor,g[,i],method='spearman')
	}
}
write.table(allcor, output, row.names = F, col.names = coln)