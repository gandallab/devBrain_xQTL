#! /usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
chr=args[1]
cutoff=as.numeric(args[2])
gname=read.table("/u/project/gandalm/cindywen/isoform_twas/gbat/out/gene_pos.txt",header=T,sep="\t")
gmap=read.table("/u/project/gandalm/cindywen/isoform_twas/gbat/data/temp_gene_mappability.txt",header=T,row.names=1)
genelist=read.table("/u/project/gandalm/cindywen/isoform_twas/gbat/out/gene_pos.txt",header=T,sep="\t")
num_reg=NULL
reg_all=NULL
mgenes=NULL
mcor=NULL
pseudo=read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/gencode.v33lift37.pseudogene.txt",as.is=T)


for(sub in 1:10){
r2=read.table(paste("/u/project/gandalm/cindywen/isoform_twas/gbat/out/pearsonR_chr",chr,"_sub",sub,".txt",sep=""),header=T,sep="\t")
cor2=read.table(paste("/u/project/gandalm/cindywen/isoform_twas/gbat/out/pearsonR_chr",chr,"_sub",sub,".txt",sep=""),header=T,sep="\t")
if(file.exists(paste("/u/project/gandalm/cindywen/isoform_twas/gbat/out/cor_chr",chr,"_sub",sub,".txt",sep=""))){
dat=read.table(paste("/u/project/gandalm/cindywen/isoform_twas/gbat/out/cor_chr",chr,"_sub",sub,".txt",sep=""),header=T,check.names=F)
r2.flag=which(is.element(r2[,1],colnames(dat)))
r2=r2[r2.flag,]
cor2=cor2[r2.flag,]
temp=as.matrix(dat)
inter.flag=which(abs(temp)<cutoff,arr.ind=T)
temp.tab=table(inter.flag[,2])
flag=as.numeric(names(which(temp.tab>0)))

r2.order=cbind(r2,1:nrow(r2))
r2.genes=merge(r2.order,gname,by.x=1,by.y=1,all.x=T)
r2.genes=r2.genes[order(r2.genes[,3]),]
m.flag=NULL
for(i in flag){
	reg.genes=gname[inter.flag[which(inter.flag[,2]==i),1],]
	target.genes=r2.genes[i,]
	if(is.na(gmap[as.character(target.genes$gene),1]) | length(which(is.element(pseudo[,1],as.character(target.genes$gene)))==T)>0) next
	if(gmap[as.character(target.genes$gene),1]>0.8){
	
		rm.flag=which(!is.element(reg.genes$gene,genelist$gene)| (reg.genes$chr==target.genes$chr & ( (reg.genes$start>target.genes$start-10000000 & reg.genes$start<target.genes$end+10000000)|(reg.genes$end>target.genes$start-10000000 & reg.genes$end<target.genes$end+10000000) )))
		rflag=NULL
		if(length(rm.flag)>0){
			if(length(rm.flag)==nrow(reg.genes)) next
			reg.genes=reg.genes[-rm.flag,]
			for(j in 1:nrow(reg.genes)){
				if(is.na(gmap[as.character(reg.genes$gene[j]),1])| length(which(is.element(pseudo[,1],as.character(reg.genes$gene[j])))==T)>0) next
				if(gmap[as.character(reg.genes$gene[j]),1]>0.8){
					rflag=c(rflag,j)
				}
				
			}
			}else{
			for(j in 1:nrow(reg.genes)){
				if(is.na(gmap[as.character(reg.genes$gene[j]),1])| length(which(is.element(pseudo[,1],as.character(reg.genes$gene[j])))==T)>0) next
				if(gmap[as.character(reg.genes$gene[j]),1]>0.8){
					rflag=c(rflag,j)
				}
				
			}
		}
		if(length(rflag)>0){
		reg_all=c(reg_all,paste(reg.genes$gene[rflag],collapse=';'))}else{
		reg_all=c(reg_all,"NA")}
		num_reg=c(num_reg,length(rflag))
		m.flag=c(m.flag,i)
	}
}
mgenes=rbind(mgenes,r2.genes[m.flag,c(1,4:6,2)])
mcor=rbind(mcor,cor2[m.flag,1:2])

}
}
chr.sig=cbind(mgenes,mcor[,2],num_reg,reg_all)
mart=read.csv("../../mart_export.txt",header=T)
mart=mart[,-c(1,2)]
temp=merge(chr.sig,mart,by.x=1,by.y=1,all.x=T)
dup.flag=which(duplicated(temp$genes))
temp=temp[-dup.flag,]
keep.flag=which(temp$num_reg>0)
print(length(keep.flag))
write.table(temp[keep.flag,],paste("chr",chr,"_sig_fdr10.txt",sep=""),row.names=F,sep="\t",quote=F)