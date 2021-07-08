library(SmartSVA)

make_sva <- function(expr,K,X){

	mod0 <- model.matrix( ~1, data = data.frame(expr) )

	if( !missing( X ) ){
		mod	<- model.matrix( ~1+X, data = data.frame(expr) )
		out	<- smartsva.cpp( dat = t(expr), mod = mod, n.sv = K, mod0 = mod0 )
	} else {
		out	<- smartsva.cpp( dat = t(expr), mod = mod0, n.sv = K )
	}

	svs	<- out$sv
	if( all( svs == 0 ) & is.null(dim(svs)) )
		svs	<- NA

	svs

}