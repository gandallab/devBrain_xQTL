#! /usr/bin/env Rscript
library(data.table)
library(argparser)
library(sva)

p <- arg_parser("")
p <- add_argument(p, "--prefix", help="Leafcutter")
p <- add_argument(p, "--geno_subj_dir", help="Subject lists")
args <- parse_args(p)

# 0: load data
dat <- fread(paste0(args$prefix, ".gz"), data.table=F)

# 1: remove duplicated header lines
dat <- dat[dat$'#Chr'!='#Chr',]

# 2: Uniform subject ID with genotype
subj <- colnames(dat)[5:ncol(dat)]
for(i in 1:length(subj)){
    subj[i]<-strsplit(subj[i],"[.]")[[1]][1]
}
index <- which(colnames(dat)=="1707.1.STARAligned.sortedByCoord.out.bam")
subj[index-4]<-"1707.1"
colnames(dat)<-c("#Chr","start","end","ID",subj)

# 3: ComBat
no_combat <- dat[,5:ncol(dat)]
rownames(no_combat) <- dat$ID
write.table(no_combat, paste0(args$prefix, "_fixSubj.tsv"), quote=F,sep="\t",col.names=T,row.names=T)

walker <- read.table(paste0(args$geno_subj_dir, "geno_subj_walker.txt"))
obrien <-read.table(paste0(args$geno_subj_dir, "geno_subj_obrien.txt"))
werling <-read.table(paste0(args$geno_subj_dir, "geno_subj_werling.txt"))
hdbr <- read.table(paste0(args$geno_subj_dir, "geno_subj_hdbr.txt"))
libd <- read.table(paste0(args$geno_subj_dir, "geno_subj_libd.txt"))

data.batch <- c()

for (i in 1:ncol(no_combat)) {
  sample <- colnames(no_combat)[i]
  if (sample %in% walker[,1]) {
    data.batch[i] <- 1
  }
  if (sample %in% obrien[,1]) {
    data.batch[i] <- 2
  }
  if (sample %in% werling[,1]) {
    data.batch[i] <- 3
  }
  if (sample %in% hdbr[,1]) {
    data.batch[i] <- 4
  }
  if (sample %in% libd[,1]) {
    data.batch[i] <- 5
  }
}

no_combatMat <- as.matrix(no_combat)
storage.mode(no_combatMat) <- "numeric"
combat <- ComBat(dat = no_combatMat, batch = data.batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
combat <- as.data.frame(combat)
write.table(combat,paste0(args$prefix, "_fixSubj_combat.tsv"), quote=F,sep="\t",col.names=T,row.names=T)

combat_bed <- cbind(dat[,1:4],combat)
write.table(combat_bed,paste0(args$prefix, "_fixSubj_combat.bed"), quote=F,sep="\t",col.names=T,row.names=F)