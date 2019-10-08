
qval <- function(xx, yy, pijmat, alpha, gamma) {
	sum1 <- 0
	sum2 <- 0

	for (yi in nrow(gamma)) {
		for (xj in ncol(gamma)) {
			val <- xx[xj]*pijmat[yi, xj]*log(alpha[yi]*gamma[yi, xj])
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

gamma_ij <- function(ii, jj, pijmat, alpha, gamma, xx, yy) {
	num <- xx[jj]*pijmat[ii, jj] + yy[ii, jj]
	xsum <- 0
	for (xj in 1:ncol(xx)) {
		xsum <- xsum + xx[xj]*pijmat[ii, xj] + yy[ii, xj]
	}
	return(num/xsum)
}

alpha_i <- function(ii, C, pijmat, alpha, gamma, xx) {
	xsum <- 0
	for (xj in 1:ncol(xx)) {
		frac <- xx[xj]/C
		elem <- frac * pijmat[ii, xj]
		xsum <- xsum + elem
	}
	return(xsum)
}

main <- function(unk=1, iters=5, converged=10e-6) {
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

	# augment y with an empty row
	emptyrow <- rep(0, nn)
	ymat <- rbind(ymat, emptyrow)
	# # aug last row w/ resids
	# resid <- (alpha %*% gamma) - beta
	# gamma_resid <- resid / alpha[kk+1]
	# gamma[kk+1,] <- gamma[kk+1,] - gamma_resid

	temp_gamma <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))
	temp_alpha <- rep(0, length(alpha))
	pijmat <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))

	print('Unknown init as:')
	print('Initial error:')
	print(sum(abs(alpha_true-alpha)))

	qhist <- rep(0, iters)

	# do em
	it <- 1
	while (it <= iters) {
		for (ii in 1:(kk+1)) {
			for (jj in 1:nn) {
				pijmat[ii, jj] <- pij(ii, jj, alpha, gamma)
			}
		}

		for (ii in 1:(kk+1)) {
			for (jj in 1:nn) {
				gij <- gamma_ij(ii, jj, pijmat, alpha, gamma, xmat, ymat)
				temp_gamma[ii, jj] = gij
			}
		}
		for (ii in 1:(kk+1)) {
			ai <- alpha_i(ii, C, pijmat, alpha, gamma, xmat)
			temp_alpha[ii] <- ai
		}
		alpha <- temp_alpha
		gamma <- temp_gamma

		qhist[it] <- qval(xmat, ymat, pijmat, alpha, gamma)
		if (it > 1 && (qhist[it] - qhist[it-1]) > 0) {
			print(paste(it, qhist[it] - qhist[it-1]))
		}
		if (it %% 100 == 0) {
			print(paste(it, '/', iters))
		}
		it <- it + 1
	}

	res <- data.frame(true=alpha_true, feast=alpha, err=abs(alpha_true-alpha))

	print(res)
	print('Total error:')
	print(sum(abs(alpha_true-alpha)))
	print('Non-unk error:')
	print(sum(abs(alpha_true[1:kk]-alpha[1:kk])))
	plot(c(1:it), qhist[1:it], 'l')
	return(list(alpha, gamma, sources, sink))
}

result <- main(iters=1000, converged=10e-6)
