
qval <- function(xx, yy, pij, alpha, gamma, clip_zero=10e-12) {
	xterm <- 0
	yterm <- 0
	for (ii in 1:nrow(gamma)) {
		for (jj in 1:ncol(gamma)) {
			logarg <- alpha[ii]*gamma[ii, jj]
			if (logarg < clip_zero) logarg <- clip_zero

			term <- xx[jj]*pij[ii, jj]*log(logarg)
			xterm <- xterm + term
			# if (is.nan(term)) print(paste(ii, jj, gamma[ii, jj]))
		}
	}

	for (ii in 1:(nrow(gamma)-1)) {
		logarg <- gamma[ii,]
		logarg[logarg < clip_zero] <- clip_zero

		term <- yy[ii,] %*% log(gamma[ii,])
		yterm <- yterm + term
	}

	return(xterm + yterm)
}

pijmat <- function(alpha, gamma) {
	pij <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))
	for (ii in 1:nrow(gamma)) {
		for (jj in 1:ncol(gamma)) {
			num <- alpha[ii]*gamma[ii, jj]
			if (num == 0) val <- 0
			else
				val <- num / (alpha %*% gamma[, jj])
			pij[ii, jj] <- val
		}
	}
	return(pij)
}

gamma_ij <- function(ii, jj, pij, alpha, gamma, xx, yy) {
	num <- xx[jj]*pij[ii, jj] + yy[ii, jj]
	denom <- t(xx) %*% pij[ii,] + sum(yy[ii,])
	term <- num/denom
	return(term)
}

alpha_i <- function(ii, C, pij, gamma, xx) {
	asum <- 1/C * (t(xx) %*% pij[ii,])
	return(asum)
}

em <- function(
	sink, sources, unk=1, iters=1000,
	converged=10e-6,
	clip_zero=10e-12,
	alpha_true=F) {

	#####################################################
	# Initialization
	#####################################################
	kk <- nrow(sources)-unk
	nn <- ncol(sources)

	ymat <- as.matrix(sources)[1:kk,]
	xmat <- as.matrix(sink)
	beta <- xmat / sum(xmat)
	C <- sum(xmat)

	alpha <- rep(1/(kk+unk), kk+unk)

	gamma <- ymat / rowSums(ymat)
	gamma[gamma < clip_zero] <- clip_zero # allocate some prob to zero entries
	gamma <- gamma / rowSums(gamma)

	if (unk > 0) {
		# Add unknown row to gamma
		# FIXME: use Liat's method
		unkrow <- rep(1/nn, nn)
		gamma <- rbind(gamma, unkrow) # recalc proportional amounts
		gamma <- gamma / rowSums(gamma)

		# Add empty row to ymat for ease of computation
		emptyrow <- rep(0, nn)
		ymat <- rbind(ymat, emptyrow)
	}

	temp_gamma <- matrix(, nrow=nrow(gamma), ncol=ncol(gamma))
	temp_alpha <- rep(0, length(alpha))


	#####################################################
	# EM
	#####################################################
	pij <- pijmat(alpha, gamma)
	print(paste('Initial Q', qval(xmat, ymat, pij, alpha, gamma)))

	it <- 1
	qhist <- c()
	while(it <= iters) {

		#####################################################
		# Compute p(i|j)
		#####################################################

		pij <- pijmat(alpha, gamma)


		#####################################################
		# Compute gamma and alpha
		#####################################################

		for (ii in 1:nrow(gamma)) {
			for (jj in 1:ncol(gamma)) {
				temp_gamma[ii, jj] <- gamma_ij(
					ii, jj,
					pij, alpha, gamma, xmat, ymat)
			}
		}
		temp_gamma[temp_gamma < clip_zero] <- clip_zero
		for (ii in 1:length(alpha)) {
			temp_alpha[ii] <- alpha_i(ii, C, pij, gamma, xmat)
		}
		gamma <- temp_gamma
		ad <- sum(abs(alpha - temp_alpha))
		alpha <- temp_alpha

		#####################################################
		# Calc Q (qnow)s and other metrics
		#  qd - difference in Q from t-1 to t
		#  r2 - r2 score between known and inferred alpha
		#  ad - total change in alpha from t-1 to t
		#####################################################

		qnow <- qval(xmat, ymat, pij, alpha, gamma)
		qhist <- c(qhist, qnow)
		qd <- 0
		if (it > 1) qd <- qhist[it] - qhist[it-1]

		r2 <- (cor(alpha, alpha_true))^2

		print(sprintf(
			'%d Q:%.2f qd:%.2f ad:%.5f, r2:%.4f',
			it, qnow, qd, ad, r2))
		# print(alpha)
		it <- it + 1

		if (ad <= converged) break
	}

	return(list(
		alpha=alpha,
		gamma=gamma,
		qhist=qhist,
		r2=r2))
}
