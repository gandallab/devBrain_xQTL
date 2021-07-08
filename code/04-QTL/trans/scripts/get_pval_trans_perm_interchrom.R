#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
chr <- as.numeric(args[1])
trans_dir <- args[2]
allp <- NULL
allrm <- NULL
gene_pos <- read.table(paste(trans_dir,"gene_pos.txt",sep = ""), header = T, as.is = T, fill = T, sep = "\t")

for(i in 1:10){
    if(file.info(paste(trans_dir,"/cor_chr",chr,"_sub",i,".txt",sep=""))$size>3){
        dat <- read.table(paste(trans_dir,"/cor_chr",chr,"_sub",i,".txt",sep=""),header=T,as.is=T)
        h2 <- read.table(paste(trans_dir,"chr",chr,"_sub",i,"_perm_h2g.txt",sep=""),header=T)
        neg.flag <- which(h2$h2g < 0 | h2$h2g > 1)
        neg.genes <- as.character(h2$gene[neg.flag])
        gene.names <- colnames(dat)
        gene.flag <- which(is.element(gene.names,neg.genes))
        if(length(gene.flag) > 0){
            dat <- data.frame(dat[,-gene.flag])
            gene.names <- gene.names[-gene.flag]
        }
        for(nc in 1:ncol(dat)){
            gene <- gene.names[nc]
            flag <- which(gene_pos$gene == gene)
            flag.chr <- gene_pos$chr[flag]

            rm.flag <- which(gene_pos$chr == flag.chr)
            allrm <- c(allrm,length(rm.flag))
            allp <- c(allp,as.numeric(dat[-rm.flag,nc]))
        }
    }
}
write.table(allp,paste(trans_dir,"/chr",chr,"_allp_h2g_filtered_interchrom.txt",sep = ""), quote = F, row.names = F, col.names = F)


i <- chr
pvals <- NULL
temp <- read.table(paste(trans_dir,"/chr",chr,"_allp_h2g_filtered_interchrom.txt",sep=""), header=T)
pvals <- c(pvals,as.numeric(temp[,1]))

observed <- sort(pvals)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

png(paste(trans_dir,"/chr",chr,"_cor_h2g_filtered_interchrom.png",sep=""))
plot(c(0,50), c(0,50), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,50), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black")
dev.off()