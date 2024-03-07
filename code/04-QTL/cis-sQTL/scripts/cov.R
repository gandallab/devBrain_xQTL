#! /usr/bin/env Rscript
library(data.table)
library(argparser)

p <- arg_parser("Generate covariates with various number of HCP")
p <- add_argument(p, "--num_hcp", help="Number of HCP in covariates")
p <- add_argument(p, "--expr", help="Normalized, ComBat")
p <- add_argument(p, "--picard", help="Compiled picard metrics")
p <- add_argument(p, "--geno_pc", help="Ancestry specific gPC")
p <- add_argument(p, "--meta", help="Metadata")
p <- add_argument(p, "--num_gpc", help="Number of gPC in covariates")
p <- add_argument(p, "--outdir", help="Out dir")

args <- parse_args(p)
num_hcp <- as.numeric(args$num_hcp)
num_gpc <- as.numeric(args$num_gpc)

# Load data
datExpr <- fread(args$expr, data.table=F)
rownames(datExpr) <- datExpr$V1
datExpr <- datExpr[,-1]

datSeq <- read.table(args$picard, header=T, stringsAsFactors=FALSE)
seq_row <- rownames(datSeq)
for (i in 1:length(seq_row)) {
  seq_row[i] <- strsplit(seq_row[i], split=" ")[[1]][1]
}
row.names(datSeq) <- seq_row
#note 1707.walker and 1707.hdbr in datSeq
rownames(datSeq)[which(rownames(datSeq) == "1707.hdbr")] <- "1707.1"
rownames(datSeq)[which(rownames(datSeq) == "1707.walker")] <- "1707"
datSeq <- datSeq[colnames(datExpr),]

geno_pc <- read.table(args$geno_pc)

meta <- read.table(args$meta, header=T, stringsAsFactors=FALSE)

# HCP
standardize<- function(X)
{
  X = as.matrix(X)
  # n = dim(X)[1]
  # p = dim(X)[2]

  X = scale(X, center = TRUE, scale = TRUE)
  # X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))

  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  return (X)
}

hidden_convariate_linear <- function(F,Y,k,lambda,lambda2,lambda3,iter) {
  ## Use Example
  # hcp = hidden_convariate_linear(standardize(datSeq), standardize(t(datExpr)),k=10,iter = 100)
  #
  # function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
  # input:
  #      F: a matrix nxd of known covariates, where n is the number of
  #      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
  #      Y: a matrix of nxg of expression data (must be standardized (columns
  #      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
  #      k: number of inferred hidden components (k is an integer)
  #      lambda, lambda2, lambda3 are model parameters
  #      (optional) iter: number of iterations (default = 100);
  #
  #      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
  #      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
  #      typically, if lambda>5, then hidden factors match the known covariates closely.
  #
  # objective:
  #
  # this function solves the following problem:
  # argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
  #
  # output:
  #      Z: matrix of hidden components, dimensionality: nxk
  #      B: matrix of effects of hidden components, dimensionality: kxg
  #      o: value of objective function on consecutive iterations.
  #
  # to use the residual data: Residual = Y - Z*B
  library(MASS)
  library(pracma)

  tol = 1e-6;

  U = matrix(0, nrow=dim(F)[2],k)
  Z = matrix(0, nrow=dim(F)[1],k)
  B = matrix(runif(dim(Z)[2]*dim(Y)[2]), nrow=dim(Z)[2], ncol=dim(Y)[2])
  F = as.matrix(F)

  n1 = dim(F)[1]
  d1 = dim(F)[2]

  n2 = dim(Y)[1]
  d2 = dim(Y)[2]

  if(n1!=n2)    stop("number of rows in F and Y must agree")

  if (k<1 | lambda<1e-6 | lambda2<1e-6 | lambda3<1e-6 ) {
    stop("lambda, lambda2, lambda3 must be positive and/or k must be an integer");
  }

  o = vector(length=iter)

  for (ii in 1:iter) {
    o[ii] = sum((Y - Z%*%B)^2) + sum((Z -  F%*%U)^2)*lambda + (sum(B^2))*lambda2 + lambda3*(sum(U^2));
    Z = (Y %*% t(B) + lambda * F %*%U) %*% ginv(B %*% t(B) + lambda * diag(dim(B)[1]))
    B = mldivide(t(Z) %*% Z + lambda2 * diag(dim(Z)[2]), (t(Z) %*% Y))
    U = mldivide(t(F) %*% F * lambda + lambda3 * diag(dim(U)[1]), lambda * t(F) %*% Z)

    if(ii > 1 &&  (abs(o[ii]-o[ii-1])/o[ii]) < tol)  break
  }

  error =  sum((Y - Z%*%B)^2) / sum(Y^2)  + sum((Z - F%*%U)^2)/sum((F%*%U)^2)
  error1 = sum((Y - Z%*%B)^2) / sum(Y^2);
  error2 = sum((Z - F%*%U)^2) / sum((F%*%U)^2);

  dz = Z%*%(B%*%t(B) + lambda*diag(dim(B)[1]))-(Y%*%t(B) + lambda*F%*%U);
  db = (t(Z)%*%Z + lambda2*diag(dim(Z)[2]))%*%B - t(Z)%*%Y;
  du = (t(F)%*%F*lambda + lambda3*diag(dim(U)[1]))%*%U-lambda*t(F)%*%Z;


  dataout = list(Z = Z, B = B, U = U)
  return(dataout)
}

# Error in svd(A) : infinite or missing values in 'x'
# scale a column with all 1s give NaN
# filter for columns with non-zero variance
datSeq <- Filter(var, datSeq)


set.seed(123)
hcp <- hidden_convariate_linear(standardize(datSeq), standardize(t(datExpr)), lambda=5,lambda2=1, lambda3=1, k=num_hcp, iter=100)
hcpMat <- cbind(1, hcp$Z) 


# Generate covariates (age + sex + gPC + HCP)
n_subj <- ncol(datExpr)
cov <- matrix(nrow = 2 + num_hcp + num_gpc , ncol = n_subj + 1)
cov <- as.data.frame(cov)

colnames(cov) <- c("id",colnames(datExpr))
for (i in 1:num_gpc) {
  cov[i,1] <- paste0("PC",i)
}
cov[num_gpc + 1, 1] <- "sex"
for (i in (num_gpc + 2):(num_gpc + 1 + num_hcp)) {
  cov[i,1] <- paste0("HCP",i-1-num_gpc)
}
cov[nrow(cov),1] <- "age"

# gPC
geno_pc <- geno_pc[,-1]
rownames(geno_pc) <- geno_pc[,1]
geno_pc <- geno_pc[,-1]
# first 2504 are 1kg
geno_pc <- geno_pc[2505:nrow(geno_pc),]

for (i in 2:ncol(cov)){
  sample <- colnames(cov)[i]
  for (j in 1:num_gpc) {
    cov[j,sample] <- geno_pc[sample,j]
    }
}

# sex and age
row.names(meta) <- meta[,1]
for (i in 2:ncol(cov)){
  sample <- colnames(cov)[i]
  cov[nrow(cov),sample] <- meta[sample,"Age"]
  cov[1 + num_gpc, sample] <- meta[sample,"inferSex"]
}

# HCP
hcp_cov <- hcpMat[,2:(num_hcp+1)]
hcp_cov <- as.data.frame(hcp_cov)
for (i in 2:ncol(cov)) {
  sample <- colnames(cov)[i]
  for(j in (num_gpc + 2):(num_gpc + 1 + num_hcp)) {
    cov[j,sample] <- hcp_cov[sample,j-1-num_gpc]
  }
}


write.table(cov, paste0(args$outdir, num_hcp, "hcp_cov.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
