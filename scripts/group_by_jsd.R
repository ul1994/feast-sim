
library('permute')
source('metrics.R')

saved <- readRDS('jsd.rda')
data <- readRDS('liat.rds')
scores <- saved[,3]

# this has some NAs...
where.min(is.na(scores))
scores[is.na(scores)] <- -1
sample_range <- 2000
targets <- c(0.125, 0.25, 0.5, 0.75)
# targets <- c(0.25)
track_index <- matrix(0, nrow=length(targets), ncol=2)
ord <- order(scores)

# assuming scores are spread nicely, find ranges to sample from
starting <- 1
for (ii in 1:length(ord)) {
	ind <- ord[ii]
	if (scores[ind] == -1) {
		starting <- starting + 1
		next
	}

	ent <- saved[ind,]
	for (ti in 1:length(targets)) {
		target <- targets[ti]
		prev <- track_index[ti, 2]
		if (abs(scores[ind] - target) < abs(prev - target)) {
			track_index[ti, 1] <- ii
			prev <- track_index[ti, 2] <- scores[ind]
		}
	}
}

print(track_index)
# print(saved[track_index[1,1]-1,])

collect_rows <- function(inds, limitN=1000) {
	mat <- matrix(, nrow=length(inds), ncol=limitN)
	for (ii in 1:length(inds)) {
		mat[ii,] <- data[inds[ii],1:limitN]
	}
	return(mat)
}

# sample around the indicies found
for (ti in 1:length(targets)) {
	target <- targets[ti]
	median <- track_index[ti, 1]
	check <- shuffle(-sample_range:sample_range)

	group <- c()
	running_avg <- -1
	for (offset in check) {
		index <- ord[median + offset]
		entry <- saved[index,]

		test_group <- group
		for (candidate in entry[1:2]) {
			candidate <- as.integer(candidate)
			if (!(candidate %in% test_group)) {
				test_group <- c(test_group, candidate)
			}
		}

		test_avg <- jsdavg(collect_rows(test_group))
		if (running_avg == -1 || abs(test_avg - running_avg) < 0.01) {
			# acceptable jsd
			# print(test_avg)
			group <- test_group
			running_avg <- test_avg
		}

		if (length(group) >= 21) break
	}

	group <- group[1:21]
	final_avg <- jsdavg(collect_rows(group))
	print(paste('Obtained:', length(group), final_avg))
	saveRDS(
		collect_rows(group),
		file=sprintf(
			'saved/jsd/sources_jsd_0%d_0%d.rds',
			as.integer(target * 1000),
			as.integer(final_avg * 1000)))
}
