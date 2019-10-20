
######################################################################
#' many_unknowns.R
#'
#' Inference with many sources
#'
#' Sample command:
#'
######################################################################

library('stats')
source('./sim.R')
source('./src.R')
source('./sims/utils.R')
source('./metrics.R')

######################################################################
#' 0. Simulation arguments
#'
#' Fixed configs
#' @param sources_file # of mixing proportions to test
#'
#' CLI parameters
#' @param T2_alphas # of mixing proportions to test
#' @param iterations Number of max EM iterations
#' @param numUnk # of unknowns to be used
#' @param maxUnk Given sources are assumed to be (K+U) many.
#'  To perform comparable tests K+1, K+2, ..., K+U the total num
#'  of U must be decided beforehand.
#'
######################################################################

# Sources file: using one fixed source
sources_file <- 'saved/unk/sources_jsd_0900_090164.rds'

args = commandArgs(trailingOnly=TRUE)

# iterations
T2_alphas <- as.integer(args[1])
print(paste('# Mixes to test:', T2_alphas))
iters <- as.integer(args[2])
print(paste('# Max iterations per test:', iters))

numUnk <- as.integer(args[3])
maxUnk <- as.integer(args[3])

######################################################################
# 1. Draw K + 2 samples S1, . . . , SK+1, SK+2 from a selected data set.
######################################################################

# Load the sources saved for having the stated JSD
raw_sources <- readRDS(sources_file)
print(paste('Total JSD is:', jsdavg(raw_sources))) # sanity check

print(sprintf('Num sources (K+%d): %d', numUnk, nrow(raw_sources))) # sanity check
print(paste('# Knowns:', (nrow(raw_sources)-maxUnk)))
print(paste('# Unknowns:', numUnk))
print(paste('# Sources given:', nrow(raw_sources)))

######################################################################
# 2. Draw noisy realization of S1, . . . , SK+1 from the
#  Multinomial distribution (denoted S^k).
######################################################################

sources <- raw_sources
# sources <- noisy_sources(raw_sources)
# print(paste('Noised JSD is:', jsdavg(sources)))

normed <- sources / rowSums(sources)
u1 <- normed[nrow(sources) - 1,]
u2 <- normed[nrow(sources),]
print(paste('JSD of last two sources:', jsd(u1, u2)))

######################################################################
# 3. For each i = 1 : T2 (different mixing proportions):
######################################################################

collected_results <- list()
for (ii in 1:T2_alphas) {

	# (a) Generate random mixing m ∼ P areto(α > 0), where Pm = 1.
	alpha_true <- generate_alphas(1,
		numK=nrow(sources)-1,
		unk=1)

	# (b) Set the sink sample abundances to m*S per taxa
	sink <- t(sources) %*% alpha_true

	# (c) Estimate the known source proportions in the sink using
	#  S^1, . . . , S^K.
	# print(paste('True unknown proportion:', alpha_true[nrow(sources)]))
	print(paste('Mix', ii, '/', T2_alphas))

	# last elements are always assumed to be the unknown sources
	#  if any provided (they are removed for FEAST)
	known_sources <- sources[1:(nrow(sources)-maxUnk),]
	results <- official_feast_wrapper(sink, known_sources, iters, unk=numUnk)

	# Keep all results for final scoring
	results$alpha_true <- alpha_true
	collected_results <- append(collected_results, list(results))

	pltfile <- sprintf('plots/unk/%d.jpg', ii)
	jpeg(pltfile, width = 350, height = 350)
	plot(c(1:length(results$qhist)), results$qhist, 'l')
	dev.off()
}

exit(0)

######################################################################
# 4. Calculate the squared Pearson correlation (r2) between the
#  estimated and the true mixing proportions per source and average
#  across sources.
######################################################################

inf_alphas <- matrix(, nrow=T2_alphas, ncol=nrow(sources))
true_alphas <- matrix(, nrow=T2_alphas, ncol=nrow(sources))
for (ti in 1:T2_alphas) {
	print(paste(
		collected_results[[ti]]$alpha[21],
		collected_results[[ti]]$alpha_true[21]))
	inf_alphas[ti,] <- collected_results[[ti]]$alpha
	true_alphas[ti,] <- collected_results[[ti]]$alpha_true
}
r2_bysource <- c()
for (si in 1:nrow(sources)) {
	r2 <- cor(inf_alphas[,si], true_alphas[,si])^2
	if (!is.na(r2)) {
		r2_bysource <- c(r2_bysource, r2)
	}
}

print(paste('For JSD:', jsdavg(sources)))
print(r2_bysource)
print(paste('Average R2:', sprintf('%.4f', mean(r2_bysource))))

######################################################################
# 5. Calculate the average Jensen-Shannon divergence of m (based on
# the pairwise Jensen-Shannon divergence).
######################################################################

# all_jsd <- c()
# for (result in collected_results) {
# 	val <- jsd(result$alpha, result$alpha_true)
# 	all_jsd <- c(all_jsd, val)
# }
# print(paste('For JSD:', jsdavg(sources)))
# print(paste('Average JSD of m:', sprintf('%.4f', mean(all_jsd))))
