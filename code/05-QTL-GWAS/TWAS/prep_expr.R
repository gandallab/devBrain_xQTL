# Generate expression file for EUR for FUSION
# Same format as GEUVADIS example data
# Regressed by covariates with 60HCP, which maximizes #eGene
# Exclude genes w/o cis-variants

# 1. Regress expression
library(data.table)
datExpr <- fread("gene.counts.scaled.normalized.bed.gz",data.table=F)
info <- datExpr[,1:4]
datExpr <- datExpr[,5:289]
cov <- read.table("60hcp_cov.txt",header=T,stringsAsFactors=F,check.names=F)

rownames(cov) <- cov$id
cov <- cov[,-1]
cov[cov=="F"] <- 0
cov[cov=="M"] <- 1
cov[cov=="Walker"] <- 1
cov[cov=="Obrien"] <- 2
cov[cov=="Werling"] <- 3
cov[cov=="HDBR"] <- 4
cov[cov=="LIBD"] <- 5
cov <- t(cov)
cov <- cbind(1, cov)
storage.mode(cov) <- "numeric"

Y <- as.matrix(datExpr)
X <- as.matrix(cov)
beta <- (solve(t(X)%*%X)%*%t(X))%*%t(Y)
datExpr.regressed <- Y - t(X[,-1]  %*% beta[-1,])

dat <- cbind(info,datExpr.regressed)
dat$end <- dat$ID
colnames(dat)[1:4]<-c("Chr","Coord","TargetID","Gene_Symbol")
dat2<-dat[,c(3,4,1,2,5:289)]

write.table(dat2,"~/project-gandalm/isoform_twas/twas/data/eur_gene_exp_regressed_60hcp.txt",col.names=T,row.names=F,quote=F,sep="\t")

# Also need an ID file
# note FID is the same as IID, not 0. Because genotype bim file is like that

# 2. Exclude genes that do not have cis-variants
eur <- fread("~/project-gandalm/isoform_twas/eqtl/results/eur.perm_60HCP/perm.all.txt.gz",data.table=F)
# > dim(eur)
# [1] 30942    11
# this number is the same for mixed, eur, amr, afr

eur <- eur[complete.cases(eur),]
# > dim(eur)
# [1] 30779    11
# this number is slightly different between ancestries, as the genotype file is different

# Load regressed expression
expr <- fread("~/project-gandalm/isoform_twas/twas/data/eur_gene_exp_regressed_60hcp.txt",data.table=F)
# > dim(expr)
# [1] 31947   289
expr2 <- expr[expr$TargetID %in% eur$V1,]
# > dim(expr2)
# [1] 30779   289
write.table(expr2,"~/project-gandalm/isoform_twas/twas/data/eur_gene_exp_regressed_60hcp.txt",col.names=T,row.names=F,quote=F,sep="\t")
