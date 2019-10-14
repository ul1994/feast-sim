
# main.R
# Main simulation study
#  R^2 score is reported for a specified JSD among sources

library('stats')
source('./sim.R')
source('./src.R')
source('./sims/utils.R')
source('./metrics.R')

######################################################################
# 0. Simulation arguments
#
# "In the simulations presented we used T1 = 10, T2 = 30, K = 20."
######################################################################

T2_alphas <- 30 # 30 mixing proportions

args = commandArgs(trailingOnly=TRUE)

# JSD
jsd_arg <- args[1]
print(paste('Choosing JSD:', jsd_arg))

# iterations
iters <- as.integer(args[2])
print(paste('# iterations per sim:', iters))


# Sources file
sources_file <- args[3]
print(paste('Loading sources:', sources_file))

######################################################################
# 1. Draw K + 1 samples S1, . . . , SK+1, from a selected data set.
######################################################################

# Load the K+1 sources saved for having such JSD
raw_sources <- readRDS(sources_file)
print(paste('Num sources (K):', nrow(raw_sources))) # sanity check
print(paste('Avg JSD is:', jsdavg(raw_sources))) # sanity check

######################################################################
# 2. Draw noisy realization of S1, . . . , SK+1 from the
#  Multinomial distribution (denoted S^k).
######################################################################

# FIXME: do this for each alpha test
sources <- noisy_sources(raw_sources)

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
	results <- em(sink, sources,
		iters=iters,
		converged=10e-6,
		alpha_true=alpha_true)
	results$alpha_true <- alpha_true

	# Keep all results for final scoring
	collected_results <- append(collected_results, list(results))

	# Save plot of Q history
	pltfile <- sprintf('plots/%s_%d.jpg', sources_file, ii)
	jpeg(pltfile, width = 350, height = 350)
	plot(c(1:length(results$qhist)), results$qhist, 'l')
	dev.off()
}

######################################################################
# 4. Calculate the squared Pearson correlation (r2) between the
#  estimated and the true mixing proportions per source and average
#  across sources.
######################################################################

all_r2 <- c()
for (result in collected_results) {
	r2score <- result$r2
	all_r2 <- c(all_r2, r2score)
}
print(paste('For JSD:', jsdavg(sources)))
print(paste('Average R2:', sprintf('%.4f', mean(all_r2))))

######################################################################
# 5. Calculate the average Jensen-Shannon divergence of m (based on
# the pairwise Jensen-Shannon divergence).
######################################################################

# print(sources)
