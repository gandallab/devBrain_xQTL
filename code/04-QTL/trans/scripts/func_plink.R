# Convert a PLINK format file (--recode 12) into a genotype matrix
# Input = f_ped,f_map strings pointing to ped/map files
# Output = W matrix, map & ped lists contining files, fam matrix containing sample data

map = read.table( f_map , stringsAsFactors=F , fill = T )
ped = read.table( f_ped , stringsAsFactors=F ) 

N = length(ped[,1])
M = length(map[,1])
W.all = matrix(0,nrow=N,ncol=M)

fam = ped[,1:6]

for ( s in 1:M ) {
    pos = c(6+s*2-1,6+s*2)
    al = unique( as.integer(c(ped[,pos[1]],ped[,pos[2]])) )
    al = al[ al != 0 ]
    for ( r in pos ) {
      W.all[ ped[,r] != al[1], s ] = W.all[ ped[,r] != al[1], s ] + 1
    }
    W.all[ ped[,r] == 0, s ] = NA
}

W = W.all
# cleanup
rm(s,W.all,al,pos,r)