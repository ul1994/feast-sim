
qval <- function(xx, yy, alpha, gamma) {
	sum1 <- 0
	sum2 <- 0

	for (yi in nrow(gamma)) {
		for (xj in ncol(gamma)) {
			val <- xx[xj]*pij(yi, xj, alpha, gamma)*log(alpha[yi]*gamma[yi, xj])
			sum1 <- sum1 + val
		}
	}

	for (yi in nrow(gamma)-1) { # FIXME: make num-unknowns changeable
		for (xj in ncol(gamma)) {
			val <- yy[yi, xj]*log(gamma[yi, xj])
			sum2 <- sum2 + val
		}
	}

	return(sum1 + sum2)
}

pij <- function(ii, jj, alpha, gamma) {
	num <- alpha[ii] * gamma[ii, jj]
	denom <- alpha %*% gamma[,jj]
	return(num/denom)
}

gamma_ij <- function(ii, jj, alpha, gamma, xx, yy) {
	yij <- if (ii > nrow(yy)) 0 else yy[ii, jj]
	num <- xx[jj]*pij(ii, jj, alpha, gamma) + yij
	xsum <- 0
	for (xj in 1:ncol(xx)) {
		yij <- if (ii > nrow(yy)) 0 else yy[ii, xj]
		# FIXME: is the +y inside or outside the sum?
		xsum <- xsum + xx[xj]*pij(ii, xj, alpha, gamma) + yij
	}
	return(num/xsum)
}

alpha_i <- function(ii, C, alpha, gamma, xx) {
	xsum <- 0
	for (xj in 1:ncol(xx)) {
		frac <- xx[xj]/C
		num <- alpha[ii]*gamma[ii, xj]
		denom <- alpha %*% gamma[, xj]

		elem <- frac * num / denom
		xsum <- xsum + elem
	}
	return(xsum)
}

main <- function(unk=1, iters=5) {
	print(paste('Init with unknowns:', unk))

	load(file="sources.Rda")
	load(file="alpha_true.Rda")
	load(file="sink.Rda")

	# prep data
	kk <- nrow(sources)-unk
	nn <- ncol(sources)
	ymat <- as.matrix(sources)[1:kk,]
	xmat <- as.matrix(sink)
	C <- sum(xmat)

	alpha <- rep(1/(kk+1), kk+1)
	gamma <- ymat / rowSums(ymat)
	beta <- xmat / sum(xmat)

	# augment gamma with a unknown source
	unkrow <- rep(1/nn, nn)
	gamma <- rbind(gamma, unkrow)
	# recalc proportional amounts
	gamma <- gamma / rowSums(gamma)

	temp_gamma <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))
	temp_alpha <- rep(0, length(alpha))
	pijmat <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))

	print('Unknown init as:')
	# print(unkrow)
	print('Initial error:')
	print(sum(abs(alpha_true-alpha)))

	qhist <- rep(0, iters)

	# do em
	for (it in 1:iters) {

		for (ii in 1:(kk+1)) {
			for (jj in 1:nn) {
				gij <- gamma_ij(ii, jj, alpha, gamma, xmat, ymat)
				temp_gamma[ii, jj] = gij
			}
		}

		gamma[1:kk,] <- temp_gamma[1:kk,]
		bhat <- alpha %*% gamma
		resid <- bhat - beta
		gamma_resid <- resid / alpha[kk+1]
		gamma[kk+1,] <- gamma[kk+1,] - gamma_resid

		for (ii in 1:(kk+1)) {
			ai <- alpha_i(ii, C, alpha, gamma, xmat)
			temp_alpha[ii] <- ai
		}
		alpha <- temp_alpha

		qhist[it] <- qval(xmat, ymat, alpha, gamma)
		if (it > 1 && qhist[it-1] > qhist[it]) {
			print(paste(it, 'Q decreased'))
			print(paste(qhist[it-1], '->', qhist[it]))
			break
		}

		if (it %% 100 == 0) {
			print(paste(it, '/', iters))
		}
	}

	res <- data.frame(true=alpha_true, feast=alpha, err=abs(alpha_true-alpha))

	print(res)
	print('Total error:')
	print(sum(abs(alpha_true-alpha)))
	print('Non-unk error:')
	print(sum(abs(alpha_true[1:kk]-alpha[1:kk])))
	plot(c(1:iters), qhist, 'l')
	return(list(alpha, gamma, sources, sink))
}

result <- main(iters=500)
