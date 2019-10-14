
loadByJSD <- function(jsdval) {
	jsdval <- sprintf('0%d', as.integer(as.numeric(jsdval) * 1000)) # stringify
	source_files <- Sys.glob(sprintf("saved/jsd/%s/sources_jsd_*.rds", jsdval))

	# just choose first in folder
	fl <- source_files[1]

	return(list(
		sources=readRDS(file=fl)
	))

	# tag <- unlist(strsplit(fl, 'sources_'))[2]
	# sources <- readRDS(file=sprintf('saved/jsd/%s/sources_%s', jsdval, tag))
	# return(list(
	# 	sources=sources
	# ))
}
