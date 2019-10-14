
# Run first:
#  find_jsd_pairs.R
#  group_by_jsd.R

library('stats')
source('./sim.R')
source('./src.R')
source('./metrics.R')

run_jsd <- '0750'
fl <- "saved/jsd/0750/sources_jsd_0750_074986.rds"


tag <- unlist(strsplit(fl, 'sources_'))[2]
alpha_true <- readRDS(file=sprintf("saved/jsd/%s/alphas_%s", run_jsd, tag))[1,]

sources <- readRDS(file=sprintf('saved/jsd/%s/sources_%s', run_jsd, tag))
sink <- readRDS(file=sprintf('saved/jsd/%s/sink_%s', run_jsd, tag))
# print(jsdavg(sources))
# sources[sources == 0] <- 1

results <- em(
	sink, sources,
	unk=0,
	iters=20, alpha_true=alpha_true)
# all_r2 <- c(all_r2, results[[4]])
# print(tag)
# print(proc.time() - ptm)