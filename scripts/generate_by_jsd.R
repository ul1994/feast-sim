

# Run first:
#  find_jsd_pairs.R
#  group_by_jsd.R

source('metrics.R')
source('sim.R')

run_jsd <- '0250'

source_files <- Sys.glob(sprintf("saved/jsd/%s/sources_*.rds", run_jsd))
for (fl in source_files) {
	tag <- unlist(strsplit(fl,'sources_'))[2]

	sources <- readRDS(fl)

	alphas <- generate_alphas(1, nrow(sources), unk=0)
	saveRDS(alphas, file=sprintf('saved/jsd/%s/alphas_%s', run_jsd, tag))

	sink <- mix_sink(alphas[1,], sources)
	saveRDS(sink, file=sprintf('saved/jsd/%s/sink_%s', run_jsd, tag))
}