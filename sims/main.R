
######################################################################
# main.R
#
# Main simulation study
#  R^2 score is reported for sources with a known JSD
#
# Sample command:
# Rscript sims/main.R 0.08 1000 saved/jsd/0080/sources_jsd_080_08108.rds
# Rscript sims/main.R 0.125 1000 saved/jsd/0125/sources_jsd_0125_012600.rds
# Rscript sims/main.R 0.5 1000 saved/jsd/0500/sources_jsd_0500_050252.rds
# Rscript sims/main.R 0.75 1000 saved/jsd/0750/sources_jsd_0750_075025.rds
# Rscript sims/main.R 0.90 1000 saved/jsd/0900/sources_jsd_0900_090009.rds
# Rscript sims/main.R (testing jsd) (iterations) (sources file to load)
######################################################################

library('stats')
source('./sim.R')
source('./src.R')
source('./sims/utils.R')
source('./metrics.R')

######################################################################
#' 0. Simulation arguments
#'
#' "In the simulations presented we used T1 = 10, T2 = 30, K = 20."
#' Fixed configs
#' @param FEAST_OFFICIAL Uses official FEAST EM if true
#' @param T2_alphas Uses official FEAST EM if true
#'
#' CLI parameters
#' @param jsd JSD currently testing for
#' @param iterations Number of max EM iterations
#' @param source_file RDS of saved sources known to have some JSD
######################################################################

FEAST_OFFICIAL <- T
T2_alphas <- 2 # 30 mixing proportions

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

# Load the K+1 sources saved for having the stated JSD
raw_sources <- readRDS(sources_file)
print(paste('Num sources (K):', nrow(raw_sources))) # sanity check
print(paste('Avg JSD is:', jsdavg(raw_sources))) # sanity check

######################################################################
# 2. Draw noisy realization of S1, . . . , SK+1 from the
#  Multinomial distribution (denoted S^k).
######################################################################

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
	print(paste('Mix', ii, '/', T2_alphas))
	if (FEAST_OFFICIAL) {
		known_only <- sources[1:(nrow(sources)-1),]
		results <- official_feast_wrapper(sink, known_only, iters)
	}
	else {
		results <- em(sink, sources,
			unk=1,
			iters=iters,
			converged=10e-6,
			alpha_true=alpha_true)

		# Save plot of Q history
		pltfile <- sprintf('plots/%s_%d.jpg', sources_file, ii)
		jpeg(pltfile, width = 350, height = 350)
		plot(c(1:length(results$qhist)), results$qhist, 'l')
		dev.off()
	}

	# Keep all results for final scoring
	results$alpha_true <- alpha_true
	collected_results <- append(collected_results, list(results))
}

######################################################################
# 4. Calculate the squared Pearson correlation (r2) between the
#  estimated and the true mixing proportions per source and average
#  across sources.
######################################################################

inf_alphas <- matrix(, nrow=T2_alphas, ncol=nrow(sources))
true_alphas <- matrix(, nrow=T2_alphas, ncol=nrow(sources))
for (ti in 1:T2_alphas) {
	inf_alphas[ti,] <- collected_results[[ti]]$alpha
	true_alphas[ti,] <- collected_results[[ti]]$alpha_true
}
r2_bysource <- c()
for (si in 1:nrow(sources)) {
	r2 <- cor(inf_alphas[,si], true_alphas[,si])^2
	r2_bysource <- c(r2_bysource, r2)
}

print(paste('For JSD:', jsdavg(sources)))
print(r2_bysource)
print(paste('Average R2:', sprintf('%.4f', mean(r2_bysource))))

# ######################################################################
# # 5. Calculate the average Jensen-Shannon divergence of m (based on
# # the pairwise Jensen-Shannon divergence).
# ######################################################################

# all_jsd <- c()
# for (result in collected_results) {
# 	val <- jsd(result$alpha, result$alpha_true)
# 	all_jsd <- c(all_jsd, val)
# }
# print(paste('For JSD:', jsdavg(sources)))
# print(paste('Average JSD of m:', sprintf('%.4f', mean(all_jsd))))
