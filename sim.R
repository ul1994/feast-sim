
library(stats)
library(actuar)

collect_sources <- function(data, numK=20, unk=1) {
	numtot <- numK + unk
	inds <- sample(1:nrow(data), numtot, replace=F)
	sources <- matrix(, nrow=numtot, ncol=ncol(data))
	it <- 1
	for (ii in inds) {
		sources[it,] <- data[ii,]
		it <- it + 1
	}
	return(sources)
}

noisy_sources <- function(counts) {
	noisy <- matrix(, nrow=nrow(counts), ncol=ncol(counts))
	for (ii in 1:nrow(counts)) {
		vector <- counts[ii,]
		vector[is.na(vector)] <- 0
		prob <- vector / sum(vector)
		sample <- rmultinom(1, size=sum(vector), prob=prob)
		noisy[ii,] <- sample
	}
	return(noisy)
}

mix_sink <- function(alpha, sources) {
	sink <- rep(0, ncol(sources))

	for (kk in 1:nrow(sources)) {
		sink <- sink + alpha[kk] * sources[kk, ]
	}

	return(sink)
}

generate_alphas <- function(batch, numK, unk=1) {
	amat <- matrix(, nrow=batch, ncol=numK+unk)
	for (ii in 1:batch) {
		alph <- rpareto(numK + unk, 3, 1)
		amat[ii,] <- alph / sum(alph)
	}
	return(amat)
}

generate_data <- function(numK=20, numAlpha=30) {
	alphamat <- generate_alphas(numAlpha, numK)
	saveRDS(alphamat, file='saved/alphas.Rda')

	raw_sources<- collect_sources(data, numK)
	for (test_i in 1:numAlpha) {
		noisy <- noisy_sources(raw_sources)
		saveRDS(noisy, file=sprintf('saved/sources%s.Rda', test_i))

		saveRDS(mix_sink(alphamat[test_i,], noisy),
			file=sprintf('saved/sink%s.Rda', test_i))
	}
}