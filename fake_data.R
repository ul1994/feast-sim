
library(stats)
library(actuar)

generate <- function(kk, nn, unk=1) {

	print(paste('N taxa (= size of sink):', nn))
	print(paste('K sources:', kk))
	print(paste('  + Unknowns:', unk))

	sources <- data.frame(
		matrix(
			vector(), kk+unk, nn,
			dimnames=list(c(), sapply(1:nn, function(tnum) paste('T', tnum, sep="")))),
		stringsAsFactors=F)

	for (val in 1:(kk+1)) {
		# sinks are characterized by sparse peaks (negative binom)
		sources[val,] <- rnbinom(nn, 20, mu=100)
	}

	# mixture is characterized by few concentrations (pareto)
	alpha_true <- rpareto(kk+unk, 3, 1)
	# manually set augment ratio
	# alpha_true[kk+1] <- sum(alpha_true[1:kk])/(1-unk_ratio)*unk_ratio
	alpha_true <- alpha_true / sum(alpha_true)

	sink <- data.frame(
		matrix(
			vector(), 1, nn,
			dimnames=list(c(), sapply(1:nn, function(tnum) paste('T', tnum, sep="")))),
		stringsAsFactors=F)
	sink <- rep(0, nn)
	# sink <- data.frame(sink=sink)

	for (val in 1:(kk+1)) {
		# simulate the sink
		sink <- sink + sources[val,] * alpha_true[val]
	}
	sink <- as.integer(sink)

	print(paste('Unknown proportion:', alpha_true[kk+1]))

	save(sources, file="sources.Rda")
	save(alpha_true, file="alpha_true.Rda")
	save(sink, file="sink.Rda")

	return(list(sources, alpha_true, sink))
}

blob <- generate(8, 128)