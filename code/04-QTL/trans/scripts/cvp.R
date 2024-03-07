source_here <- function(x, ...) {
    dir <- "."
    if(sys.nframe()>0) {
        frame <- sys.frame(1)
        if (!is.null(frame$ofile)) {
            dir <- dirname(frame$ofile)
        }
    }
    source(file.path(dir, x), ...)
}
source_here("func_reml.R")

# uncomment "solve" function below to override default matrix inversion
#  -- if standard solve fails to calculate invese of input, calculate a pseudoinverse
# solve <- base::solve  # to go back to normal
## solve <- function(M, epsilon=1e-8){
##     standardMatrixInverseWorked <- F
##     inverseM <- NA
##     try({
##         inverseM <-  base::solve(M)
##         standardMatrixInverseWorked <- T
##     }, silent=T)
##     if(!standardMatrixInverseWorked){
##         udvt <- svd(M)
##         d <-  (udvt$d)
##         di <- rep(0,length(d))
##         di[abs(d)>epsilon] <- 1/d[abs(d)>epsilon]
##         inverseM <- udvt$v %*% diag(di) %*% t(udvt$u)
##     }
##     return(inverseM)   
## }

ones <- function(n,p){ # n rows p columns of ones
    #ret <- rep(1,n) %*% t(rep(1,p))
    ret <- matrix(1,n,p)
    return(ret)
}

tl <- function(m, n=5, p=5){
    ret <- m[1:n, 1:p]
    return(ret)
}

setTermWidth <- function(){ options(width=Sys.getenv("COLUMNS")) }


 
## getBlups <- function(y, X=NA, Z, sig2g=NA, sig2e=NA, h2=NA, return_bBlupsLoo=F, return_intermediates=F)
## alternative cv scheme. fixed effects estimated (loocv),
## then given those, residuals (y-Xbeta_loocv) are calculated
## gBLUPs then calculated using loocv blup formula and residuals
## (intercept is updated and cross-validated with cvBLUPs )
## calculate cvBLUPs and related quantities.
## y = Xbeta + Zb + epsilon,
## cvblups,blups from linear mixed model with 1 structured variance component + noise
## give h2 or (sig2g and sig2e)
## X = data.frame or matrix of covariates for fixed effects (not including intercept), intercept added in function
## y = phenotype/outcome
## Z = data.frame or matrix of (scaled, centered) genotype matrix: n rows, m columns for n subjects, m snps
##
## if(return_bBlupsLoo=T): return bBlupsLoo
## bBlupsLoo := m by n matrix of bBLUPs (BLUP effect estimates) corresponding to n by m genotype matrix (but transposed)
## column j,  bBlupsLoo[,j] is a set of m bBlups calculated with subject j removed
## row i, bBlupsLoo[i,] gives effect estimates corresponding to column i of scaled genotype matrix Z
##
## if(return_intermediates=T): return intermediates := list(   V ,  Vi ,  H ,  sig2g ,  sig2e   )
## fixed effects loocv estimated once, then given the fixed effects, the gBLUPs are estimated by loocv
## sig2g,sig2e are estimated once using all data,
## then when each subject is left out in turn,  bBlups are reestimated and gBlupLoo predicitons for left-out subjects made
## betaFe is the GLS estimate of beta, yfe is X %*% betaFe, yfeLoo is the loocv estimate of Xbeta
## bBlup is the vector of standard BLUP effect size estimates
## gBlup are the standard in-sample BLUP fits
## gBlupLoo are the loocv out of sample gBLUP predictions
## pc is the prior probability of a snp being causal ( >0, <=1, 1 for infinitessimal model, <1 spike-slab)
getBlups <- function(y, X=NA, Z, sig2g=NA, sig2e=NA, h2=NA, return_bBlupsLoo=F, return_intermediates=F){
    y <- as.numeric(y)  # funny problems if y is a scale() return object
    Z <- as.matrix(Z) # not data.frame
    n <- nrow(Z)
    m <- ncol(Z)  
    K <- 1/m * tcrossprod(Z)

    if(all(is.na(X))){
        X <- as.matrix(data.frame(rep(1,n))) # force to be a matrix, not vector
    }else{                                            # R is wierd
        X <- as.matrix(data.frame(cbind(rep(1,n),X))) # force to be a matrix, not vector, not wide, not data.frame
    }
 
    if(is.na(sig2g) || is.na(sig2e)){
        betaOLS <- solve(crossprod(X)) %*% (t(X) %*% y)
        XbetaOLS <- as.numeric(X %*% betaOLS)
        sig2g.temp <- (var(y-XbetaOLS))*h2
        sig2e.temp <- (var(y-XbetaOLS))*(1-h2)
        Vinv.temp <- solve(sig2g.temp*K + sig2e.temp*as.matrix(diag(n))  )
        betaGLS <-   (solve(t(X) %*% Vinv.temp %*% X) %*% (t(X) %*% Vinv.temp %*% y))
        XbetaGLS <- as.numeric(X %*% betaGLS)
        sig2g <- (var(y-XbetaGLS))*h2
        sig2e <- (var(y-XbetaGLS))*(1-h2)
    }
    
 
    V <- sig2g*K + sig2e*as.matrix(diag(n))
    #save(list = ls(), file="inCVBLUPfunction.rdata")
    Vi <- solve(V)

    xtvixixtvi <- solve(t(X) %*% Vi %*% X) %*% (t(X) %*% Vi)
    betaFe <- xtvixixtvi %*% y 
    Xproj <- X %*% xtvixixtvi
    yfe <- Xproj %*% y
    yfeLoo <- (yfe - diag(Xproj)*y)/(1-diag(Xproj))
    #P <- as.matrix(diag(n)) - Xproj # projection operator to move to orthog space of X
    yTilde <- y - yfeLoo # phenotype residuals, fixed effect subtracted
    yFeResid <- y - yfe # not cross-validated fixed effect

    # fix the intercept, recenter for each fold in cross-validation
    X1 <- matrix(1,n,1)
    X1proj <- X1 %*%solve(t(X1) %*% Vi %*% X1) %*% (t(X1) %*% Vi)
    P1 <-  as.matrix(diag(n)) - X1proj
    
    bBlupStandard <- (1/m) * sig2g * t(Z) %*% (Vi %*% yFeResid)
    bBlup <- (1/m) * sig2g * t(Z) %*% (Vi %*% P1 %*% yTilde)
    H <- sig2g * K %*% Vi %*% P1 
    gBlup <- H %*% yTilde
    gBlupLoo <- (gBlup - diag(H)*yTilde) / (1 - diag(H) )

    ## if(pc < 1){
    ##     bShrink <- doShrink(bBlups=bBlup, m=m, pc=pc, n=n, sig2g=sig2g, sig2e=sig2e){
    ## }
    
    ret <- list(
               betaFe=betaFe, yfe=yfe,  yfeLoo=yfeLoo,
               bBlup=bBlup, bBlupStandard=bBlupStandard,
               gBlup=gBlup, gBlupLoo=gBlupLoo
    )

    if(return_bBlupsLoo){
        bBlupsLoo <- getBBlupsLoo(y=y,  Z=Z, Vi=Vi , yfeLoo=yfeLoo , gBlupLoo=gBlupLoo, sig2g=sig2g  )
        ret$bBlupsLoo <- bBlupsLoo
    }

    if(return_intermediates){
        intermediates <- list(   V=V,  Vi=Vi,  H=H,  sig2g=sig2g,  sig2e=sig2e  )
        ret$intermediates <- intermediates
    }
    
   
    return(ret)
       
}


test_getBlups <- function(seed=1,  useX=T){
    set.seed(seed)
    sig2g <- 0.5
    sig2e <- 0.5
    sig2fe <- 0.2
    n <- 1000
    m <- 1000
    p <- 5  # number of fixed effect covariates
    pc <- 1  # fraction of snps that are causal
    mc <- ceiling(pc*m)
    mnc <- m - mc
    

    beta0 <- rnorm(1)  # intercept
    beta <- rnorm(p) * sqrt(sig2fe/p)
    b <- c( rnorm(mc) * sqrt(sig2g/mc) , rep(0, mnc))
    
    X <-  matrix(rnorm(n*p),n,p) 
    
    G <- matrix(rbinom(n*m,size=2,prob=.5),n,m)
    
    Z <- as.matrix(scale(G))
    
    e <- rnorm(n)*sqrt(sig2e)
     
    Xbeta <- as.numeric(X %*% beta)
    
    Zb <- as.numeric(Z %*% b)
     
    y <- beta0 + Xbeta + Zb + e
    

    if(!useX){
        X <- NA
    }
    
    blups <- getBlups(y=y, X=X, Z=Z, sig2g=sig2g, sig2e=sig2e)
    blupsOld <- getBlupsOld(y=y, X=X, Z=Z, sig2g=sig2g, sig2e=sig2e)
    bblups <- getBlups(y=y, X=X, Z=Z, sig2g=sig2g, sig2e=sig2e, return_bBlupsLoo=T)
    blupsAndIntermediates <- getBlups(y=y, X=X, Z=Z, sig2g=sig2g, sig2e=sig2e, return_bBlupsLoo=T)
    ret <-  list(blups=blups, blupsOld=blupsOld, bblups=bblups, blupsAndIntermediates=blupsAndIntermediates)

    return(ret)  
}

 
## getFixed <- function(y, X=NA, sig2g=NA, sig2e=NA, h2=NA)
## y = Xbeta + Zb + epsilon
## X  fixed effect covariates n-by-p
## a column of ones will be added to X for intercept term
## If X is NA, will just use an intercept term
## Z scaled genotypes n-by-m
## If K is not given, build it given Z
## beta fixed effects, b ~N(0, sig2g/m Im), epsilon ~ N(0,sig2e In)
## betaHat is the GLS estimate of beta, yfe is X %*% betaFe,
## yfeLoo is the loocv estimate of Xbeta
## give sig2g and sig2e or h2.  h2 will be converted to sig2g and sig2e (using var(y))
getFixed <- function(y, X=NA, Z=NA, K=NA, sig2g=NA, sig2e=NA, h2=NA){
    y <- as.numeric(y)  # funny problems if y is a scale() return object

    
    if(all(is.na(K))){
        Z <- as.matrix(Z) # not data.frame
        m <- ncol(Z)
        n <- nrow(Z)
        K <- 1/m * tcrossprod(Z)
    }else{
        K <- as.matrix(K)
        n <- nrow(K)
    }

    
    if(all(is.na(X))){
        X <- as.matrix(data.frame(rep(1,n))) # force to be a matrix, not vector
    }else{                                            # R is weird
        X <- as.matrix(data.frame(cbind(rep(1,n),X))) # force to be a matrix, not vector, not wide, not data.frame
    }
 
    if(is.na(sig2g) || is.na(sig2e)){
        betaOLS <- solve(crossprod(X)) %*% (t(X) %*% y)
        XbetaOLS <- as.numeric(X %*% betaOLS)
        sig2g.temp <- (var(y-XbetaOLS))*h2
        sig2e.temp <- (var(y-XbetaOLS))*(1-h2)
        Vinv.temp <- solve(sig2g.temp*K + sig2e.temp*as.matrix(diag(n))  )
        betaGLS <-   (solve(t(X) %*% Vinv.temp %*% X) %*% (t(X) %*% Vinv.temp %*% y)) # basic GLS calculation 
        XbetaGLS <- as.numeric(X %*% betaGLS)
        sig2g <- (var(y-XbetaGLS))*h2
        sig2e <- (var(y-XbetaGLS))*(1-h2)
    }
    
 
    V <- sig2g*K + sig2e*as.matrix(diag(n))
    Vi <- solve(V)

    xtvixixtvi <- solve(t(X) %*% Vi %*% X) %*% (t(X) %*% Vi)
    betaFe <- xtvixixtvi %*% y 
    Xproj <- X %*% xtvixixtvi
    yfe <- Xproj %*% y
    yfeLoo <- (yfe - diag(Xproj)*y)/(1-diag(Xproj))
     
    
    return(list( betaFe=betaFe, yfe=yfe,  yfeLoo=yfeLoo ) )
}

test_getFixed <- function(seed=1){
    set.seed(seed)

    sig2g <- 0.5
    sig2e <- 0.5
    sig2fe <- 0.2
    n <- 500
    m <- 500
    p <- 5  # number of fixed effect covariates
    pc <- 1  # fraction of snps that are causal
    mc <- ceiling(pc*m)
    mnc <- m - mc
     

    beta0 <- rnorm(1)  # intercept
    beta <- rnorm(p) * sqrt(sig2fe/p)
    b <- c( rnorm(mc) * sqrt(sig2g/mc) , rep(0, mnc))
    
    X <-  matrix(rnorm(n*p),n,p) 
     
    G <- matrix(rbinom(n*m,size=2,prob=.5),n,m)
     
    Z <- as.matrix(scale(G))
    K <- (1/m) * tcrossprod(Z)
    
    e <- rnorm(n)*sqrt(sig2e)
    
    Xbeta <- as.numeric(X %*% beta)
     
    Zb <- as.numeric(Z %*% b)
     
    y <- beta0 + Xbeta + Zb + e
     
    fe1 <- getFixed(y=y, X=X, Z=Z, sig2g=sig2g, sig2e=sig2e )
    fe2 <- getFixed(y=y, X=X, K=K, sig2g=sig2g, sig2e=sig2e )
    fe3 <- getFixed(y=y, X=NA, Z=Z, sig2g=sig2g, sig2e=sig2e )
    fe4 <- getFixed(y=y, X=X, Z=Z, h2=.5 )

    ret <- list(fe1=fe1, fe2=fe2, fe3=fe3, fe4=fe4)

    for(i in 1:4){
        for (j in (1)){
            print(var(ret[[i]][[j]]))
        }
    }

    return(ret)
}

 

# out of sample BLUP predictions
# y = Xbeta + Zb + epsilon
# given fixed and random effect estimates, predict outcomes (y) for new subjects
# number of new subjects = nnew
# yHatNew = (fixed effect prediction) + (random effect prediction)
# yHatNew = yHatNewFe + yHatNewRe
# for Z = genotypes, yHatNewRe will be the genetic predictions for new subjects
# beta:= fixed effects for covariates X
# b:= random effects for covariates (scaled genotypes) Z
# betaHat := vector of extimates of beta, length=p+1
# bBlups := vector of BLUPs for b, length=m
# Xnew := matrix of fixed effect covariates for new subjects, nnew-by-p (not including column of ones)
# Znew := matrix of random effect genotypes for new subjects, nnew-by-m
# let first column of Xnew be all 1s, and first element of betaHat be intercept term
# if Xnew or Znew are centered and scaled, center and scale using mean and SD of X,Z (reference data for betaHat, bBlups)
oosBlupPredictions <- function(bBlups,betaHat, Xnew=NA, Znew){
    n <- nrow(Znew)
    if(all(is.na(Xnew))){
        Xnew <- as.matrix(data.frame(rep(1,n))) # force to be a matrix, not vector
    }else{                                            # R is wierd
        Xnew <- as.matrix(data.frame(cbind(rep(1,n),Xnew))) # force to be a matrix, not vector, not wide, not data.frame
    }   
    yHatNewFe <- Xnew %*% betaHat
    yHatNewRe <- Znew %*% bBlups
    yHatNew <- yHatNewFe + yHatNewRe
    ret <- list(yHatNewFe=yHatNewFe, yHatNewRe=yHatNewRe, yHatNew=yHatNew)
    return(ret)
}

# get cross-validated BLUPs for n subjects given genotyeps G and covariates X
# make genetic predictions for nnew additional subjecs using the trained on the original n subjects
# y = outcomes for original n subjects
# X, Xnew fixed covariates (without column of 1s)
# G, Gnew = unscaled genotytpes for fixed effects (0,1,2; additive coding)
# returns list(gBlupLoo, gNewPred)
# gBlupLoo = cross-validated genetic predictions for original n subjects
# gNewPred = genetic predictions for new subjects
#            with genotypes Gnew and covariates Xnew based on model from (y,X,Z)
# extra=T for more stuff returned
getBlupsAndOutOfSamplePredictions <- function(y, X, G, Xnew, Gnew, extra=F, blupsFunc=getBlups){
    m <- ncol(G)
    n <- nrow(G)
    nnew <- nrow(Gnew)
    Gscale <- scale(G)
    Z <- as.matrix(Gscale)
    K <- tcrossprod(Z)/m
    Znew <- (Gnew -  ones(nnew,m)%*%diag(as.numeric(attr(Gscale, "scaled:center")))) /
        (ones(nnew,m)%*%diag(as.numeric(attr(Gscale, "scaled:scale"))))
    if(all(is.na(X))){
        Xuse <- matrix(1,n,1)
    }
    else{
        Xuse <- cbind(rep(1,n),X)
    }
    varComp <- aiREML(A=list(K), y=y, X=Xuse, Var=c(.5,.5),F) # add col of 1s to X
    #cvp <- getBlups(y=y, X=X, Z=Z, h2=varComp[[1]])
    cvp <- blupsFunc(y=y, X=X, Z=Z, h2=varComp[[1]])
    oosPredictions <- oosBlupPredictions(bBlups=cvp$bBlup, betaHat=cvp$betaFe, Xnew=Xnew, Znew=Znew)
    if(extra){
        ret <- list(gBlupLoo=as.numeric(cvp$gBlupLoo),gBlup=as.numeric(cvp$gBlup),
                yfe=as.numeric(cvp$yfe), yfeLoo=as.numeric(cvp$yfeLoo),
                gNewPred=as.numeric(oosPredictions$yHatNewRe),
                    varComp=varComp)
    }else{
        ret <- list(gBlupLoo=as.numeric(cvp$gBlupLoo), 
                    gNewPred=as.numeric(oosPredictions$yHatNewRe) )
    }
                    
    return(ret)
}

test_getBlupsAndOutOfSamplePredictions <- function(seed=1, blupsFunc=getBlups, useX=T){
    set.seed(seed)
    sig2g <- 0.5
    sig2e <- 0.5
    sig2fe <- 0.2
    n <- 1000
    m <- 1000
    p <- 5  # number of fixed effect covariates
    pc <- 1  # fraction of snps that are causal
    mc <- ceiling(pc*m)
    mnc <- m - mc
    nnew <- 100  # number of additional hold-out subjects to predict

    beta0 <- rnorm(1)  # intercept
    beta <- rnorm(p) * sqrt(sig2fe/p)
    b <- c( rnorm(mc) * sqrt(sig2g/mc) , rep(0, mnc))
    
    X <-  matrix(rnorm(n*p),n,p) 
    Xnew <- matrix(rnorm(nnew*p),nnew,p) 
    G <- matrix(rbinom(n*m,size=2,prob=.5),n,m)
    Gnew <- matrix(rbinom(nnew*m,size=2,prob=.5),nnew,m)
    Z <- as.matrix(scale(G))
    Znew <- as.matrix(scale(Gnew))
    e <- rnorm(n)*sqrt(sig2e)
    enew <- rnorm(nnew)*sqrt(sig2e)
    Xbeta <- as.numeric(X %*% beta)
    Xnewbeta <- as.numeric(Xnew %*% beta)
    Zb <- as.numeric(Z %*% b)
    Znewb <- as.numeric(Znew %*% b)
    y <- beta0 + Xbeta + Zb + e
    ynew <- beta0 + Xnewbeta + Znewb + enew

    if(!useX){
        X <- NA
        Xnew <- NA
    }
    
    cvp.oos <- getBlupsAndOutOfSamplePredictions(y=y, X=X, G=G,
                                                 Xnew=Xnew, Gnew=Gnew, extra=T, blupsFunc=blupsFunc)
    ret <- list(y=y, Zb=Zb, cvBlups=cvp.oos$gBlupLoo, ynew=ynew,
                Znewb=Znewb, gNewPred=cvp.oos$gNewPred, varComp=cvp.oos$varComp)

    cat("\n\n\ny = Xbeta + Zb + e\n")
    cat("gBlups is in-sample estimate of Zb (random effect contribution to phenotype)\n")
    cat("cvBlups is loocv estimate of Zb  \n")
    cat("yfe is in-sample estimate of Xbeta (fixed effect contribution to phenotype)\n")
    cat("yfeLoo is loocv estimate of Xbeta  \n")
    cat("gNewPred is prediciton into new subjects using betaHat and bBlups from linear mixed model\n\n")
    cat("Variance component estimation\n")
    print(cvp.oos$varComp[1:3])
    cat("\nCorrelations\n")
    print(round(cor(data.frame(y=y, Zb=Zb, gBlups=cvp.oos$gBlup, cvBlups=cvp.oos$gBlupLoo,
                               Xbeta=Xbeta,yfe=cvp.oos$yfe, yfeLoo=cvp.oos$yfeLoo,  e=e)),3))
    cat("\nhold-out sample predictions (new data)\n")
    print(round(cor(data.frame(ynew=ynew, Znewb=Znewb,
                               gNewPred=cvp.oos$gNewPred, Xnewbeta=Xnewbeta, enew=enew)),3))

    return(ret)  
}


 

# read a binary GCTA format grm with name mygrm.grm.bin from file and convert to a symmetric R matrix
# grm <- ReadGRMBin(prefix="mygrm")$GRM
# from library gap
ReadGRMBin <- function (prefix, AllN = FALSE, size = 4) 
{
    BinFileName <- paste(prefix, ".grm.bin", sep = "")
    NFileName <- paste(prefix, ".grm.N.bin", sep = "")
    IDFileName <- paste(prefix, ".grm.id", sep = "")
    id <- read.table(IDFileName)
    n <- dim(id)[1]
    BinFile <- file(BinFileName, "rb")
    grm <- readBin(BinFile, n = n * (n + 1)/2, what = numeric(0), 
        size = size)
    close(BinFile)
    NFile <- file(NFileName, "rb")
    if (AllN) 
        N <- readBin(NFile, n = n * (n + 1)/2, what = numeric(0), 
            size = size)
    else N <- readBin(NFile, n = 1, what = numeric(0), size = size)
    close(NFile)
    i <- sapply(1:n, function(i) i * (i + 1)/2)
    GRM <- matrix(NA, n, n)
    GRM[upper.tri(GRM, diag = TRUE)] <- grm
    GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
    invisible(list(grm = grm, id = id, N = N, GRM = GRM))
}

# write a symmetric matrix "K" with row/col ids "id" to file in GCTA binary GRM format with name "prefix"
writeMat2Grm <- function(prefix, K, id, N=NA){
    grm <- grm[upper.tri(K,diag=T)]
    if(!(all(is.na(N)))){
        AllN <- T
    }else{
        AllN <- F
    }
    WriteGRMBin(prefix=prefix, grm=grm, id=id, N=N, AllN=AllN)
}

# from library gap
 WriteGRMBin<- function (prefix, grm, id, size = 4, N=NA, AllN=F) 
{
    BinFileName <- paste(prefix, ".grm.bin", sep = "")
    NFileName <- paste(prefix, ".grm.N.bin", sep = "")
    IDFileName <- paste(prefix, ".grm.id", sep = "")
    grm.bin <- file(BinFileName, "wb")
    writeBin(grm, grm.bin, size = size)
    close(grm.bin)
    if(AllN){
        grm.N.bin <- file(NFileName, "wb")
        writeBin(N, grm.N.bin, size = size)
        close(grm.N.bin)
    }
    write.table(id, IDFileName, col.names = FALSE, quote = FALSE, 
        row.names = FALSE, sep = "\t")
}



## getBBlupsLoo <- function(y,  Z, Vi , yfeLoo , gBlupLoo, sig2g=NA,  h2=NA  )
##  use this if have y,Z,Vi,yfeLoo,gBlupLoo and (either sig2g or h2)
##  otherwise just use getBlups() function with getBlups( return_bBlupsLoo=T )
##
## get bBlupsLoo := m by n matrix of bBLUPs (BLUP effect estimates) corresponding to n by m genotype matrix (but transposed)
## column j,  bBlupsLoo[,j] is a set of m bBlups calculated with subject j removed
## row i, bBlupsLoo[i,] gives effect estimtes corresponding to column i of scaled genotype matrix Z
##
## y = beta0 + Xbeta + Zb + epsilon;
## beta0 mean,
## X n-by-p covarites (do not include column of 1s for intercept)
## Z n-by-m scaled genoyptes (n subjects, m SNPs)
## beta0, beta fixed effects; b random effects
## var(b) = sig2g/m;  sig2g genetic variance;  var(epsilon)= sig2e
## h2=sig2g/(sig2g+sig2e);  (sig2g + sig2e) = var(y - Xbeta)
## K = ZZt/m;  V = (sig2g*K + sig2e*I_n)
## yfeLoo = LOO cross-validated estimtes of fixed effects (Xbeta) as from getBlups() or getFixed functions
## gBlupLoo = LOO cross-validated estiamtes of genetic effect (Zb) as from getBlups() function
##
## bBlupsLoo <- getBBlupsLoo(y=y, Z=Z, sig2g=sig2g, Vi=Vi, yFeLoo=yFeLoo, gBlupLoo=gBlupLoo)
## bBlupsLoo <- getBBlupsLoo(y=y, Z=Z, h2=h2, Vi=Vi, yFeLoo=yFeLoo, gBlupLoo=gBlupLoo)  
getBBlupsLoo <- function(y,  Z, Vi , yfeLoo , gBlupLoo, sig2g=NA,  h2=NA  ){
    y <- as.numeric(y)  # funny problems if y is a scale() return object
    Z <- as.matrix(Z) # not data.frame
    n <- nrow(Z)
    m <- ncol(Z)

    yTilde <- as.numeric(y - yfeLoo)
    
    if( is.na(sig2g)  ){
        sig2g <- var(yTilde) * h2      
    }  

    Y <- as.matrix( diag(yTilde) ) %*% matrix( 1, n, n )
    diag( Y ) <- gBlupLoo
    bBlupsLoo <- sig2g/m * t(Z) %*% Vi %*% Y

    return(bBlupsLoo)
} 


#  univariate nonuniform shrink -- convert bBLUPs to predictions for Bayesian sparse linear mixed model
doShrink <- function(bBlups, pc, n, m, sig2g, sig2e, isBeta=F){
    # bBlups (vector of blups on effect size scale)
    # m features,
    # pc = prob(feature causal = non-zero effect size),
    # n observations
    # is.Beta (independent marginal GWAS betas) (do uniform shrink and posterior-causal shrink)
    # if !is.Beta (then BLUP) only do posterior-causal schrink (but adust scale for m vs pc*m in uniform shrink)
    
    u <- bBlups # old notation, u for bblups
    sig2p <- sig2g + sig2e
    h2 <- sig2g/sig2p
    d <- (sig2g/m/pc + sig2p/n)
    numi <- pc/sqrt(d) * exp(-u*u/2/d)
    deni <- numi + (1-pc)/sqrt(sig2p/n) * exp(-u*u/2*n/sig2p)
    pc.i <- numi/deni
    if(isBeta){
        shrunku <- (h2/(h2+m*pc/n)) * pc.i * u
    }
    else{
        shrunku <- pc.i * u *  (h2/(h2+m*pc/n))/ (h2/(h2+m/n))
    }
    return(shrunku)
}


## get sparse cross-validated predictions
## getSparse_CVPs <- function(pc, bBlupsLoo, gBlupsLoo, Z, y, sig2g, sig2e, returnSparseBBlups=F) 
## pc is fraction causal or a vector of k fractions causal to consider: pc=c(0.5, 0.1, 0.02), pc=.01, 
## returns a n-by-(k+2) matrix of genetic predictions (unshrunk, shrunk@pc, combinedCVP)
## the "combined" are the loocv predicitons from a linear model: y ~ (unshrunk) + (shrunk@pc[1]) + (shrunk@pc[2]) ...
## bBlupsLoo is a matrix with different sets of bBlups in each column, m-by-n for LOOCV.
## But generally could be any matrix of bBlups (each column is processed separately)
getSparse_CVPs <- function(pc, bBlupsLoo, gBlupLoo, Z, y, sig2g, sig2e, returnSparseBBlups=F){
    Z <- as.matrix( Z )
   
    m <- ncol(Z)
    n <- nrow(Z)
    k <- length(pc)

    cvpSparse <- matrix(NA, n, k+2)
    colnames(cvpSparse) <- c("gBlupLoo", paste("gBlupLoo.",pc[1:k],sep=""), "combinedCVP")
    cvpSparse[,1] <- gBlupLoo
   
    sparseBBlups <- list()

    for(j in 1:k){
        bBlupsLooSparse <- as.matrix( bBlupsLoo )
        for( i in 1:ncol(bBlupsLoo)){
            bBlupsLooSparse[,i] <- doShrink( bBlups=bBlupsLoo[,i], m=m, pc=pc[j], n=n, sig2g=sig2g, sig2e=sig2e)
        }
        cvpSparse[,j+1] <- apply(( t(bBlupsLooSparse) * Z ), 1, sum, na.rm=T)
        if( returnSparseBBlups ){
            sparseBBlups[length(sparseBBlups)+1] <- bBlupsLooSparse
            names(sparseBBlups)[length(sparseBBlups)] <- paste("bBlups.",pc[j],sep="")
        }
        
    }

    xtemp <- as.matrix(data.frame(intercept=rep(1,length(y)),cvpSparse[,1:(k+1)]))
    Hx <- xtemp %*% solve(crossprod(xtemp)) %*% t(xtemp)
    cblup <- Hx %*% y  # combined BLUP
    cblupLoo <- (cblup - diag(Hx)*y)/(1-diag(Hx))

    cvpSparse[,k+2] <- cblupLoo
    ret <- list( cvpSparse = cvpSparse  )
    if( returnSparseBBlups ){
        ret$sparseBBlups=sparseBBlups
    }
    return(ret)
       
}

test_getSparse_CVPs <- function(pcGenerate=0.02, pcFit=c(.5, .1, .01), seed=1){
    set.seed(seed)
    sig2g <- 0.5
    sig2e <- 0.5
    sig2fe <- 0.2
    n <- 1000
    m <- 1500
    p <- 5  # number of fixed effect covariates
    pc <- pcGenerate  # fraction of snps that are causal
    mc <- ceiling(pc*m)
    mnc <- m - mc
    

    beta0 <- 0 #rnorm(1)  # intercept
    beta <- rnorm(p) * sqrt(sig2fe/p)
    b <- c( rnorm(mc) * sqrt(sig2g/mc) , rep(0, mnc))
    X <-  matrix(rnorm(n*p),n,p)   
    G <- matrix(rbinom(n*m,size=2,prob=.5),n,m)
    Z <- as.matrix(scale(G))
    e <- rnorm(n)*sqrt(sig2e)
    Xbeta <- as.numeric(X %*% beta)
    Zb <- as.numeric(Z %*% b)
    y <- beta0 + Xbeta + Zb + e
    K <- tcrossprod(Z)/m
    
    varComp <- aiREML(A=list(K), y=y, X=cbind(rep(1,n),X), Var=c(.5,.5),F) # add col of 1s to X
    sig2g <- varComp$vars[1,1]
    sig2e <- varComp$vars[2,1]
    V <- sig2g*K + sig2e*diag(n)
    Vi <- solve(V)
    
    blupSet <- getBlups(y=y, X=X, Z=Z, sig2g=sig2g, sig2e=sig2e, return_bBlupsLoo=T, return_intermediates=T ) 

     
    sparse_CVPs <- getSparse_CVPs(pc=pcFit,  bBlupsLoo=blupSet$bBlupsLoo, gBlupLoo=blupSet$gBlupLoo,
                                      Z=Z, y=y, sig2g=blupSet$intermediates$sig2g, sig2e=blupSet$intermediates$sig2e, F)

    cvpSparse <- sparse_CVPs$cvpSparse
    ret <- data.frame(y=y, Zb=Zb, Xbeta=Xbeta, e=e, XbetaHatLoo=blupSet$yfeLoo,
                      cvBlups=blupSet$gBlupLoo, cvpSparse=cvpSparse)
    yfeLoo <- as.numeric(blupSet$yfeLoo)
    cat("\nCorrelations:\n")
    temp <- cbind(y,Xbeta,Zb,e, yfeLoo,cvpSparse)
    print( round( cor( temp), 3 ) )
    ## print(round(cor(a),3) )
    ## par(mfrow=c(2,2))
    ## plot(sparse_cvBlups ~ cvBlups, data=ret, pch=20, col=rgb(0,0,1,.6), xlim=c(-1.5,1.5), ylim=c(-1.5,1.5) ); abline(0,1)
    ## plot(sparse_cvBlups ~ Zb, data=ret, pch=20, col=rgb(0,0,1,.6), xlim=c(-1.5,1.5), ylim=c(-1.5,1.5) ); abline(0,1)
    ## plot( cvBlups ~ Zb, data=ret, pch=20, col=rgb(0,0,1,.6), xlim=c(-1.5,1.5), ylim=c(-1.5,1.5) ); abline(0,1)
    ## plot(sparse_cvBlups ~ y, data=ret, pch=20, col=rgb(0,0,1,.6), xlim=c(-1.5,1.5), ylim=c(-1.5,1.5) ); abline(0,1)
    
    ## abline(0,1)
    
    return(ret)  
}


  
## ## getBlupsOld <- function(y, X=NA, Z, sig2g=NA, sig2e=NA, h2=NA)
## ## calculate cvBLUPs and related quantities.
## ## y = Xbeta + Zb + epsilon,
## ## cvblups,blups from linear mixed model with 1 structured variance component + noise
## ## give h2 or (sig2g and sig2e)
## ## X = data.frame or matrix of covariates for fixed effects (not including intercept), intercept added
## ## y = phenotype/outcome
## ## Z = data.frame or matrix of (scaled, centered) genotype matrix: n rows, m columns for n subjects, m snps
## ## 
## ## fixed and random effects are cross-validated together --
## ## sig2g,sig2e are estimated once using all data,
## ## then when each subject is left out in turn, beta and bBlups are reestimated
## ## betaFe is the GLS estimate of beta, yfe is X %*% betaFe, yfeLoo is the loocv estimate of Xbeta
## ## bBlup is the vector of standard BLUP effect size estimates
## ## gBlup are the standard in-sample BLUP fits
## ## gBlupLoo are the loocv out of sample gBLUP predictions
getBlupsOld <- function(y, X=NA, Z, sig2g=NA, sig2e=NA, h2=NA){
    y <- as.numeric(y)  # funny problems if y is a scale() return object
    Z <- as.matrix(Z) # not data.frame
    n <- nrow(Z)
    m <- ncol(Z)  
    K <- 1/m * tcrossprod(Z)

    if(all(is.na(X))){
        X <- as.matrix(data.frame(rep(1,n))) # force to be a matrix, not vector
    }else{                                            # R is wierd
        X <- as.matrix(data.frame(cbind(rep(1,n),X))) # force to be a matrix, not vector, not wide, not data.frame
    }
 
    if(is.na(sig2g) || is.na(sig2e)){
        betaOLS <- solve(crossprod(X)) %*% (t(X) %*% y)
        XbetaOLS <- as.numeric(X %*% betaOLS)
        sig2g.temp <- (var(y-XbetaOLS))*h2
        sig2e.temp <- (var(y-XbetaOLS))*(1-h2)
        Vinv.temp <- solve(sig2g.temp*K + sig2e.temp*as.matrix(diag(n))  )
        betaGLS <-   (solve(t(X) %*% Vinv.temp %*% X) %*% (t(X) %*% Vinv.temp %*% y))
        XbetaGLS <- as.numeric(X %*% betaGLS)
        sig2g <- (var(y-XbetaGLS))*h2
        sig2e <- (var(y-XbetaGLS))*(1-h2)
    }
    
 
    V <- sig2g*K + sig2e*as.matrix(diag(n))
    Vi <- solve(V)

    xtvixixtvi <- solve(t(X) %*% Vi %*% X) %*% (t(X) %*% Vi)
    betaFe <- xtvixixtvi %*% y 
    Xproj <- X %*% xtvixixtvi
    yfe <- Xproj %*% y
    yfeLoo <- (yfe - diag(Xproj)*y)/(1-diag(Xproj))
    P <- as.matrix(diag(n)) - Xproj # projection operator to move to orthog space of X
    
    bBlup <- (1/m) * sig2g * t(Z) %*% (Vi %*% (P  %*% y))
    H <- sig2g * K %*% Vi %*% P
    gBlup <- H %*% y
    gBlupLoo <- (gBlup - diag(H)*y) / (1 - diag(H) )
    
    return(list(
        betaFe=betaFe, yfe=yfe,  yfeLoo=yfeLoo,
        bBlup=bBlup, gBlup=gBlup, gBlupLoo=gBlupLoo) )
}