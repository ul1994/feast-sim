
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

gamma_ij <- function(ii, jj, pij, alpha, gamma, xx, yy) {
	num <- xx[jj]*pij[ii, jj] + yy[ii, jj]
	denom <- t(xx) %*% pij[ii,] + sum(yy[ii,])
	return(num/denom)
}

alpha_i <- function(ii, C, pij, gamma, xx) {
	asum <- 1/C * (t(xx) %*% pij[ii,])
	return(asum)
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
				temp_gamma[ii, jj] <- gamma_ij(ii, jj, pij, alpha, gamma, xmat, ymat)
			}
		}
		for (ii in 1:length(alpha)) {
			temp_alpha[ii] <- alpha_i(ii, C, pij, gamma, xmat)
		}
		gamma <- temp_gamma
		alpha <- temp_alpha

		qnow <- qval(xmat, ymat, pij, alpha, gamma)
		qhist[it] <- qnow
		if (iters > 50 && it %% 100 == 0) {
			print(paste(it, 'Q:', qnow))
		}
		if (iters <= 50) {
			print(paste(it, 'Q:', qnow))
		}
		it <- it + 1
	}

	if (iters > 50) {
		res <- data.frame(true=alpha_true, feast=alpha, err=abs(alpha_true-alpha))
		print(res)
		print('Total error:')
		print(sum(abs(alpha_true-alpha)))
		print('Non-unk error:')
		print(sum(abs(alpha_true[1:kk]-alpha[1:kk])))
		plot(c(1:it), qhist[1:it], 'l')
	}

	return(list(alpha, gamma, sources, sink))
}

result <- main(iters=1000, converged=10e-6)
