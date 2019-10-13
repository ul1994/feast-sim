
qval <- function(xx, yy, pij, alpha, gamma) {
	xterm <- 0
	yterm <- 0
	for (ii in 1:nrow(gamma)) {
		for (jj in 1:ncol(gamma)) {
			term <- xx[jj]*pij[ii, jj]*log(alpha[ii]*gamma[ii, jj])
			xterm <- xterm + term
		}
	}

	for (ii in 1:(nrow(gamma)-1)) {
		term <- yy[ii,] %*% log(gamma[ii,])
		yterm <- yterm + term
	}

	return(xterm + yterm)
}

pijmat <- function(alpha, gamma) {
	pij <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))
	for (ii in 1:nrow(gamma)) {
		for (jj in 1:ncol(gamma)) {
			val <- alpha[ii]*gamma[ii, jj] / (alpha %*% gamma[, jj])
			pij[ii, jj] <- val
		}
	}
	return(pij)
}

gamma_ij <- function(ii, jj, pij, alpha, gamma, xx, yy, clip_zero=10e-3) {
	num <- xx[jj]*pij[ii, jj] + yy[ii, jj]
	denom <- t(xx) %*% pij[ii,] + sum(yy[ii,])
	term <- num/denom
	if (term < clip_zero) term <- clip_zero
	return(term)
}

alpha_i <- function(ii, C, pij, gamma, xx) {
	asum <- 1/C * (t(xx) %*% pij[ii,])
	return(asum)
}

em <- function(sink, sources, unk=1, iters=1000,
	converged=10e-6,
	clip_zero=10e-3) {

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

	# augment ymat with an empty row for ease of computation
	emptyrow <- rep(0, nn)
	ymat <- rbind(ymat, emptyrow)

	temp_gamma <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))
	temp_alpha <- rep(0, length(alpha))

	qhist <- rep(0, iters)
	it <- 1
	print('Unknown init as:')
	print('Initial Qval:')
	pij <- pijmat(alpha, gamma)
	print(qval(xmat, ymat, pij, alpha, gamma))

	while(it <= iters) {
		pij <- pijmat(alpha, gamma)

		for (ii in 1:nrow(gamma)) {
			for (jj in 1:ncol(gamma)) {
				temp_gamma[ii, jj] <- gamma_ij(
					ii, jj,
					pij, alpha, gamma, xmat, ymat)
			}
		}
		for (ii in 1:length(alpha)) {
			temp_alpha[ii] <- alpha_i(ii, C, pij, gamma, xmat)
		}
		gamma <- temp_gamma
		alpha <- temp_alpha

		qnow <- qval(xmat, ymat, pij, alpha, gamma)
		qhist[it] <- qnow
		qd <- 0
		if (it > 1) qd <- qhist[it] - qhist[it-1]
		print(sprintf(
			'%d Q:%.2f qd:%.2f',
			it, qnow, qd))

		it <- it + 1
	}

	return(list(alpha, gamma, qhist))
}
