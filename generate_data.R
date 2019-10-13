
source('sim.R')

numAlpha <- 30
numK <- 20
tiny <- T

alphamat <- generate_alphas(numAlpha, numK)
# saveRDS(alphamat, file='saved/alphas.Rda')

# data <- readRDS('liat.rds')

for (test_i in 1:numAlpha) {
	raw_sources<- collect_sources(data, numK)
	if (tiny) {
		raw_sources<- collect_sources_small(data, numK)
	}
	noisy <- noisy_sources(raw_sources)
	# saveRDS(noisy, file=sprintf('saved/sources%s.Rda', test_i))

	sink <- mix_sink(alphamat[test_i,], noisy)
	# saveRDS(,
	# 	file=sprintf('saved/sink%s.Rda', test_i))
}