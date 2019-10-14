
# Run first:
#  find_jsd_pairs.R
#  group_by_jsd.R

library('stats')
source('./sim.R')
source('./src.R')
source('./metrics.R')

unk = 1
print(paste('Init with unknowns:', unk))

source('metrics.R')
source('sim.R')

run_jsd <- '0750'
source_files <- Sys.glob(sprintf("saved/jsd/%s/sources_jsd_*.rds", run_jsd))
print(paste('Running', length(source_files), 'files...'))

all_r2 <- c()
for (fl in source_files) {
	print(fl)
	tag <- unlist(strsplit(fl, 'sources_'))[2]

	alphas <- readRDS(file=sprintf("saved/jsd/%s/alphas_%s", run_jsd, tag))[1,]
	# ptm <- proc.time()

	alpha_true <- alphas
	sources <- readRDS(file=sprintf('saved/jsd/%s/sources_%s', run_jsd, tag))
	sink <- readRDS(file=sprintf('saved/jsd/%s/sink_%s', run_jsd, tag))
	print(jsdavg(sources))

	# results <- em(sink, sources, iters=40, alpha_true=alpha_true)
	# all_r2 <- c(all_r2, results[[4]])
	# print(tag)
	# print(proc.time() - ptm)
}
