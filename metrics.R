
jsdavg <- function(srcmat) {
	proportions <- srcmat / rowSums(srcmat)
	proportions <- srcmat / rowSums(srcmat)

	mat <- jsdmatrix(proportions)

	jsum <- c()
	for (ii in 1:(ncol(mat)-1)) {
		for (jj in (ii+1):ncol(mat)) {
			jsum <- c(jsum, mat[ii, jj])
		}
	}
	return(sum(jsum) / length(jsum))
}

jsdmatrix <- function(x){
	d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
	for(i in 1:(nrow(x)-1)){
		for(j in (i+1):nrow(x)){
			d[i,j] <- jsd(x[i,], x[j,])
			d[j,i] <- d[i,j]
		}
	}
	return(d)
}

jsd <- function(p,q){
	m <- (p + q)/2
	return((kld(p,m) + kld(q,m))/2)
}

kld <- function(p,q){
	nonzero <- p>0 & q>0
	return(sum(p[nonzero] * log2(p[nonzero]/q[nonzero])))
}
