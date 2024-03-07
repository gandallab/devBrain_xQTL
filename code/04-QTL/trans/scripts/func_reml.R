# Library required for deltamethod computation of SE
library('msm')

# Utility function for calculating log of determinant
`logdet` <-
function(p) {
	det_l <- determinant(p,log=T)

	if ( det_l$sign[1] == 1 ) return(det_l$modulus[1])
	else return(det_l$modulus[1] * -1);
}

# Performs checks that REML input is consistent
`validate_input` <-
function( A , y , Var , X = NULL ) {
	N = length(y)
	
	# Check VC lengths
	if ( length(A) + 1 != length(Var) ) {
		cat("ERROR: Number of variance components doesn't match number of items in Kinship list\n")
		return(FALSE)
	}

	# Test each matrix
	for ( i in 1:length(A) ) {
		if ( !isSymmetric(A[[i]]) ) {
			cat("ERROR: Matrix ",i," is not symmetric\n",sep='')
			if ( isSymmetric( round(A[[i]],4)) ) {
				cat("Matrix ",i," is approximately symmetric, continuing\n",sep='')
			} else {
				return(FALSE)
			}
		}
		if ( nrow(A[[i]]) != N || ncol(A[[i]]) != N ) {
			cat("ERROR: Matrix ",i," is not the correct size - ",dim(A[[i]]),"\n",sep='')
			return(FALSE)
		}
	}
	
	# Test fixed effects
	if ( !is.null(X) ) {
		if ( nrow(X) != N ) { 
			cat("ERROR: Matrix of fixed effects is incorrect size - ",dim(X),"\n",sep='')
			return(FALSE)
		}
	}
	
	return( TRUE )
}

# Execute Average-Information ML (no fixed effects) until convergence
# A = list of GRM
# y = vector of phenotype entries
# Var = initial variance components (percent)

`aiML` <-
function( A, y , Var , verbose=TRUE , CALC_SE=TRUE , BOUNDED=FALSE ){

	if ( !validate_input( A , y , Var ) ) { return(NA) }
	
	r <- length(A) + 1
	N <- length(y)

	# Add matrix of residuals to A
	A[[r]] <- diag(N)

	AI <- matrix(0, ncol=r, nrow=r)
	S <- matrix(0, ncol=r, nrow=r)
	s <- matrix(0, ncol=1, nrow=r)

	l_dif <- 10
	it <- 0

	Var <- var(y) * Var

	# Perform a single iteration of EM-based REML to initiate parameters
	if(verbose) cat("Performing a single iteration of EM-REML\n")
	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv
	logL <- -0.5 * ( logdet(V) + t(y) %*% P %*% y )

	if(verbose)
	cat("Prior values from EM-REML:")
	for ( i in 1:r ) {
		Var[i] <- (Var[i]^2 * t(y) %*% P %*% A[[i]] %*% P %*% y + sum(diag(Var[i]*diag(N) - Var[i]^2 * P %*% A[[i]])) )/N
		if(verbose)
		cat(" ",Var[i],sep='')
	}
	if(verbose)
	cat('\n')

	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv
	logL <- -0.5 * ( logdet(V) + t(y) %*% P %*% y )

	if(verbose)
	cat ("EM:\t",logL,'\n',sep='')

	# Iterate AI REML until convergence
	# while ( abs(l_dif) >= 10^-4 & it < 100 ){

	# ** GCTA style:
	while ( it < 100 & ( abs(l_dif) >= 10^-4 | (abs(l_dif) < 10^-2 & l_dif < 0)) ){

		it <- it + 1

		# Average information matrix
		for ( i in 1:r ) {
			for( ii in 1:r ) {
				if ( i == r && ii == r ) AI[r,r] <- t(y) %*% P %*% P %*% P %*% y
				else if ( i == r ) AI[r,ii] <- t(y) %*% P %*% P %*% A[[ii]] %*% P %*% y
				else if ( ii == r ) AI[i,r] <- t(y) %*% P %*% A[[i]] %*% P %*% P %*% y
				else AI[i,ii] <- t(y) %*% P %*% A[[i]] %*% P %*% A[[ii]] %*% P %*% y
			}
		}
		AI <- 0.5*AI

		# Vector of first derivatives of log likelihood  function
		for ( i in 1:r ) {
			if ( i == r ) s[r,1] <- sum(diag(( P ))) - ( t(y) %*% P %*% P %*% y )
			else s[i,1] <- sum(diag(( P %*% A[[i]] ))) - ( t(y) %*% P %*% A[[i]] %*% P %*% y )
		}
		s <- -0.5*s

		# New variance components from AI and likelihood
		# Var <- Var + solve(AI) %*% s

		# ** GCTA style:
		if ( l_dif > 1 ) Var <- Var + 0.316*(solve(AI) %*% s)
		else Var <- Var + solve(AI) %*% s

		# Re-calculate V and P matrix
		V <- 0
		for ( i in 1:r ) V <- V + A[[i]] * Var[i]
		Vinv <- solve(V)
		P <- Vinv 

		# Likelihood of the MLM (ignoring constants)
		new_logL <- -0.5 * ( logdet(V) + t(y) %*% P %*% y )
		l_dif <- new_logL - logL
		logL <- new_logL

		if(verbose) {
		cat(it,'\t',logL,sep='')
		for( i in 1:r ) cat( '\t',Var[i],sep='' )
		cat('\n')
		}

		if(BOUNDED) {
			if( min(Var/sum(Var)) < 0 ) {
				if(verbose) cat("VC has escaped parameter space, bounding and exiting\n")
				break()
			}
		}
	}
	if(verbose)
	cat('\n')

	if ( !CALC_SE ) {
		return( list( "h2" = Var[1]/sum(Var) , "vc" = Var ))

	} else {

	# Calculate matrix for standard errors (same as AI matrix but w/out y)
	for( i in 1:r ) {
		for ( ii in 1:r ) {
			S[i,ii] <- sum(diag(P %*% A[[i]] %*% P %*% A[[ii]] ))
		}
	}
	S <- 0.5*S
	Sinv <- solve(S)

	if(verbose){
	for( i in 1:r ) cat( "V(G",i,")\t",Var[i],'\t',sqrt(Sinv[i,i]),'\n',sep="")
	}
	
	# Construct string equation for delta method "~x1+x2 ... +xn"
	sum_all.eq = ""
	for(i in 1:r) {
		if ( i == 1 ) sum_all.eq <- paste(sum_all.eq,"x",i,sep='')
		else sum_all.eq <- paste(sum_all.eq,"+x",i,sep='')
	}
	SE.p <- deltamethod( as.formula(paste("~",sum_all.eq)),Var,Sinv,ses=T)

	if( verbose ) {
	cat( "Vp\t",sum(Var),'\t',SE.p,'\n',sep="")
	}

	SE.i <- rep(0,r)
	
	for( i in 1:r ) {
		# Construct string equation for delta method "~xi/(x1 ... +xn)"
		SE.eq <- paste("~x",i,"/(",sum_all.eq,")",sep='')
		SE.i[i] <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)
		
		if(verbose)
		cat( "V(G",i,")/Vp\t",Var[i]/sum(Var),'\t',SE.i[i],'\n',sep='')
	}
	
	return( list( "h2" = Var[1]/sum(Var) , "se" = SE.i , "vc" = Var ))
	}
}

# Execute binary search for ML h^2
# A = list of GRM
# y = vector of phenotype entries
# Var = initial variance components (percent)

`binsML` <-
function( A, y , Var , verbose = TRUE ){

	if ( !validate_input( A , y , Var ) ) { return(NA) }

	r <- length(A) + 1
	N <- length(y)

	# Add matrix of residuals to A
	A[[r]] <- diag(N)

	h2o = 0.5
	besth2o = h2o
	step = 0.25

	# Compute log-likelihood
	Var = c( var(y) * h2o , var(y) * (1-h2o) )
	V = 0
	for ( i in 1:r ) V = V + A[[i]] * Var[i]
	Vinv = solve(V)
	logL = -0.5 * ( logdet(V) + t(y) %*% Vinv %*% y )
	pre_bestlogl = 0
	bestlogl = logL

	it = 0
	while ( it < 15 && (abs(bestlogl - pre_bestlogl) > 10^-4 || (bestlogl == pre_bestlogl && abs(logL - pre_bestlogl) > 10^-4)) ) {
		it = it + 1
		pre_bestlogl = logL
		if(verbose) cat(it,besth2o,bestlogl,'\n')

		# Move down
		h2o = besth2o - step; 
		# Compute log-likelihood
		Var = c( var(y) * h2o , var(y) * (1-h2o) )
		V = 0
		for ( i in 1:r ) V = V + A[[i]] * Var[i]
		Vinv = solve(V)
		logL = -0.5 * ( logdet(V) + t(y) %*% Vinv %*% y )

		if(logL > bestlogl) {
			besth2o = h2o
			bestlogl = logL
			step = step / 2.0
			next()
		}

		# Move up		
		h2o = besth2o + step;
		# Compute log-likelihood
		Var = c( var(y) * h2o , var(y) * (1-h2o) )
		V = 0
		for ( i in 1:r ) V = V + A[[i]] * Var[i]
		Vinv = solve(V)
		logL = -0.5 * ( logdet(V) + t(y) %*% Vinv %*% y )
		
		if(logL > bestlogl) {
			besth2o=h2o
			bestlogl=logL
			step = step / 2.0
			next()
		}
		step = step / 2.0;
	}
	if(verbose) cat(besth2o,bestlogl,'\n')

	return( list( "h2" = besth2o ))
}

# Execute Average-Information REML until convergence
# A = list of GRM
# y = vector of phenotype entries
# X = matrix of fixed effects
# Var = initial variance components (percent)

`aiREML` <-
function( A, y, X , Var , verbose = TRUE , CALC_SE = TRUE ){

	if ( !validate_input( A , y , Var , X ) ) { return(NA) }

	r <- length(A) + 1
	N <- length(y)

	# Add matrix of residuals to A
	A[[r]] <- diag(N)

	AI <- matrix(0, ncol=r, nrow=r)
	S <- matrix(0, ncol=r, nrow=r)
	s <- matrix(0, ncol=1, nrow=r)

	l_dif <- 10
	it <- 0

	Var <- var(y) * Var

	# Perform a single iteration of EM-based REML to initiate parameters
	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv - Vinv %*% X %*% solve( t(X) %*% Vinv %*% X ) %*% t(X) %*% Vinv

	if(verbose)
	cat("Prior values from EM-REML:")
	for ( i in 1:r ) {
		Var[i] <- (Var[i]^2 * t(y) %*% P %*% A[[i]] %*% P %*% y + sum(diag(Var[i]*diag(N) - Var[i]^2 * P %*% A[[i]])) )/N
		if(verbose)
		cat(" ",Var[i],sep='')
	}
	if(verbose)
	cat('\n')

	V <- 0
	for ( i in 1:r ) V <- V + A[[i]] * Var[i]
	Vinv <- solve(V)
	P <- Vinv - Vinv %*% X %*% solve( t(X) %*% Vinv %*% X ) %*% t(X) %*% Vinv
	logL <- -0.5 * ( logdet(V) + logdet(t(X) %*% Vinv %*% X) + t(y) %*% P %*% y )

	# Calculate using fixed effects
	# b <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
	# logL <- -0.5 * ( t(y - X * b) %*% Vinv %*% (y - X * b) + logdet(V) + logdet(t(X) %*% Vinv %*% X) )

	if(verbose)
	cat ("EM:\t",logL,'\n',sep='')

	# Iterate AI REML until convergence
	# while ( abs(l_dif) >= 10^-4 & it < 100 ){

	# ** GCTA style:
	while ( it < 100 & ( abs(l_dif) >= 10^-4 | (abs(l_dif)<10^-2 & l_dif < 0)) ){

		it <- it + 1

		# Average information matrix
		for ( i in 1:r ) {
			for( ii in 1:r ) {
				if ( i == r && ii == r ) AI[r,r] <- t(y) %*% P %*% P %*% P %*% y
				else if ( i == r ) AI[r,ii] <- t(y) %*% P %*% P %*% A[[ii]] %*% P %*% y
				else if ( ii == r ) AI[i,r] <- t(y) %*% P %*% A[[i]] %*% P %*% P %*% y
				else AI[i,ii] <- t(y) %*% P %*% A[[i]] %*% P %*% A[[ii]] %*% P %*% y
			}
		}
		AI <- 0.5*AI

		# Vector of first derivatives of log likelihood  function
		for ( i in 1:r ) {
			if ( i == r ) s[r,1] <- sum(diag(( P ))) - ( t(y) %*% P %*% P %*% y )
			else s[i,1] <- sum(diag(( P %*% A[[i]] ))) - ( t(y) %*% P %*% A[[i]] %*% P %*% y )
		}
		s <- -0.5*s

		# New variance components from AI and likelihood
		# Var <- Var + solve(AI) %*% s

		# ** GCTA style:
		if ( l_dif > 1 ) Var <- Var + 0.316*(solve(AI) %*% s)
		else Var <- Var + solve(AI) %*% s

		# Re-calculate V and P matrix
		V <- 0
		for ( i in 1:r ) V <- V + A[[i]] * Var[i]
		Vinv <- solve(V)
		P <- Vinv - Vinv %*% X %*% solve( t(X) %*% Vinv %*% X ) %*% t(X) %*% Vinv

		# Likelihood of the MLM (ignoring constants)
		new_logL <- -0.5 * ( logdet(V) + logdet(t(X) %*% Vinv %*% X) + t(y) %*% P %*% y )
		l_dif <- new_logL - logL
		logL <- new_logL

		if(verbose) {
		cat(it,'\t',logL,sep='')
		for( i in 1:r ) cat( '\t',Var[i],sep='' )
		cat('\n')
		}
	}
	if(verbose)
	cat('\n')

	if ( !CALC_SE ) {
		return( list( "h2" = Var[1]/sum(Var) , "vc" = Var ))

	} else {
	# Calculate matrix for standard errors (same as AI matrix but w/out y)
	for( i in 1:r ) {
		for ( ii in 1:r ) {
			S[i,ii] <- sum(diag(P %*% A[[i]] %*% P %*% A[[ii]] ))
		}
	}
	S <- 0.5*S
	Sinv <- solve(S)

	if(verbose){
	for( i in 1:r ) cat( "V(G",i,")\t",Var[i],'\t',sqrt(Sinv[i,i]),'\n',sep="")
	}
	
	# Construct string equation for delta method "~x1+x2 ... +xn"
	SE.eq <- "~"
	for(i in 1:r) {
		if ( i == 1 ) SE.eq <- paste(SE.eq,"x",i,sep='')
		else SE.eq <- paste(SE.eq,"+x",i,sep='')
	}
	SE.p <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)
	if(verbose) {
	cat( "Vp\t",sum(Var),'\t',SE.p,'\n',sep="")
	}
	
	SE.i <- rep(0,r)

	for( i in 1:r ) {
		# Construct string equation for delta method "~xi/(x2 ... +xn)"
		SE.eq <- paste("~x",i,"/(",sep='')
		for( ii in setdiff(1:r,i) ) {
			if ( ii == setdiff(1:r,i)[1] ) SE.eq <- paste (SE.eq,"x",i,sep='')
			SE.eq <- paste (SE.eq,"+x",ii,sep='')
		}
		SE.eq <- paste(SE.eq,")",sep='')

		SE.i[i] <- deltamethod(as.formula(SE.eq),Var,Sinv,ses=T)
		if(verbose){
		cat( "V(G",i,")/Vp\t",Var[i]/sum(Var),'\t',SE.i[i],'\n',sep='')
		}
	}

	return( list( "h2" = Var[1]/sum(Var) , "se" = SE.i , "vc" = Var ))
	}
}