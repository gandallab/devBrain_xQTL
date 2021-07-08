# Standardize a genotype matrix
# Input = W matrix of genotypes
# Output = normalized matrix W, mafs vector listing frequencies, miss vector listing missingness

N = nrow( W )
M = ncol( W )

mafs = colMeans(W,na.rm=T)/2
miss = apply( W,2, function(x) sum(is.na(x)) )

# Center the markers
W.new = matrix(nrow=N,ncol=M)
for ( i in 1:N ) {
	W.new[i,] = W[i,] - 2*mafs
}

# Set missing values to zero
W.new[ is.na(W.new) ] = 0

# Re-scale the variances to 1
vars = apply(W.new,2,var)
for ( i in 1:N ) {
	W.new[i,] = W.new[i,]/sqrt(vars)
}

# Zero out monomorphics
W.new[ is.na(W.new) ] = 0
W = W.new
rm(W.new)